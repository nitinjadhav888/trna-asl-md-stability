"""
Phase 6: Backbone Root-Mean-Square Deviation (RMSD) Analysis.
Calculates backbone RMSD independently for all 6 replicates across:
  1. Full 17-nt ASL
  2. 5-bp Stem
  3. 7-nt Loop (Stem-Aligned / Global and Loop-Aligned / Internal)
Outputs:
  - analysis/tables/rmsd_replicates.csv
  - analysis/tables/rmsd_summary.csv
  - analysis/tables/rmsd_time_series.csv
"""
import os
import sys
import numpy as np
import pandas as pd
import mdtraj as md

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_DIR = os.path.join(BASE_DIR, 'analysis', 'tables')
os.makedirs(TABLES_DIR, exist_ok=True)

REPLICATES = [
    ('wt', 1), ('wt', 2), ('wt', 3),
    ('cua', 1), ('cua', 2), ('cua', 3)
]

STEM_RESIDUES = [0, 1, 2, 3, 4, 12, 13, 14, 15, 16] # residues 27-31 and 39-43
LOOP_RESIDUES = [5, 6, 7, 8, 9, 10, 11]              # residues 32-38
BACKBONE_ATOMS = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]

def run_rmsd_analysis():
    print("=" * 80)
    print("PHASE 6: BACKBONE RMSD ANALYSIS")
    print("=" * 80)
    
    replicate_rows = []
    time_series_dfs = []
    
    for construct, rep_num in REPLICATES:
        rep_tag = f"{construct}_rep{rep_num}"
        dcd_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_imaged.dcd')
        pdb_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_ref.pdb')
        
        # Load canonical starting reference from Phase 4C final equilibrated structure
        p4c_pdb = os.path.join(BASE_DIR, 'simulations', construct, 'equilibration', 'final_equilibrated.pdb')
        p4c_full = md.load(p4c_pdb)
        solute_idx = p4c_full.topology.select('not water and not (resname K or resname CL)')
        ref = p4c_full.atom_slice(solute_idx)
        
        traj = md.load(dcd_path, top=pdb_path)
        n_frames = traj.n_frames
        times_ns = np.linspace(0.010, 20.000, n_frames)
        
        # Identify atom indices for backbone subsets
        bb_all = [a.index for a in traj.topology.atoms if a.name in BACKBONE_ATOMS]
        stem_bb = [i for i in bb_all if traj.topology.atom(i).residue.index in STEM_RESIDUES]
        loop_bb = [i for i in bb_all if traj.topology.atom(i).residue.index in LOOP_RESIDUES]
        
        # 1. Full ASL RMSD (align on full backbone, evaluate full backbone)
        traj_full = traj[:]
        traj_full.superpose(ref, atom_indices=bb_all)
        diff_full = traj_full.xyz[:, bb_all, :] - ref.xyz[0, bb_all, :]
        rmsd_full = np.sqrt(np.mean(np.sum(diff_full**2, axis=-1), axis=-1)) * 10.0 # Angstroms
        
        # 2. Stem RMSD & Global Loop RMSD (align strictly on stem backbone)
        traj_stem = traj[:]
        traj_stem.superpose(ref, atom_indices=stem_bb)
        diff_stem = traj_stem.xyz[:, stem_bb, :] - ref.xyz[0, stem_bb, :]
        rmsd_stem = np.sqrt(np.mean(np.sum(diff_stem**2, axis=-1), axis=-1)) * 10.0 # Angstroms
        
        # Global loop RMSD (evaluated while aligned to stem backbone to capture hinge tilt)
        diff_loop_glob = traj_stem.xyz[:, loop_bb, :] - ref.xyz[0, loop_bb, :]
        rmsd_loop_global = np.sqrt(np.mean(np.sum(diff_loop_glob**2, axis=-1), axis=-1)) * 10.0 # Angstroms
        
        # 3. Internal Loop RMSD (align strictly on loop backbone, evaluate on loop backbone)
        traj_loop = traj[:]
        traj_loop.superpose(ref, atom_indices=loop_bb)
        diff_loop_int = traj_loop.xyz[:, loop_bb, :] - ref.xyz[0, loop_bb, :]
        rmsd_loop_internal = np.sqrt(np.mean(np.sum(diff_loop_int**2, axis=-1), axis=-1)) * 10.0 # Angstroms
        
        # Final window (last 10 ns: frames 1000 to 2000)
        final_slice = slice(1000, 2000)
        
        # Store time-series
        df_ts = pd.DataFrame({
            'frame': np.arange(n_frames),
            'time_ns': np.round(times_ns, 4),
            'construct': construct,
            'replicate': rep_num,
            'replicate_tag': rep_tag,
            'rmsd_full_asl_angstrom': np.round(rmsd_full, 4),
            'rmsd_stem_angstrom': np.round(rmsd_stem, 4),
            'rmsd_loop_global_angstrom': np.round(rmsd_loop_global, 4),
            'rmsd_loop_internal_angstrom': np.round(rmsd_loop_internal, 4)
        })
        time_series_dfs.append(df_ts)
        
        # Replicate summary statistics
        def get_stats(arr):
            return {
                'mean': float(np.mean(arr)),
                'std': float(np.std(arr)),
                'median': float(np.median(arr)),
                'p25': float(np.percentile(arr, 25)),
                'p75': float(np.percentile(arr, 75)),
                'p95': float(np.percentile(arr, 95)),
                'final_10ns_mean': float(np.mean(arr[final_slice])),
                'final_10ns_std': float(np.std(arr[final_slice]))
            }
            
        st_full = get_stats(rmsd_full)
        st_stem = get_stats(rmsd_stem)
        st_loop_glob = get_stats(rmsd_loop_global)
        st_loop_int = get_stats(rmsd_loop_internal)
        
        row = {
            'construct': construct,
            'replicate': rep_num,
            'replicate_tag': rep_tag,
            'frames': n_frames,
            'duration_ns': 20.0,
            # Full ASL
            'full_asl_mean': st_full['mean'],
            'full_asl_std': st_full['std'],
            'full_asl_median': st_full['median'],
            'full_asl_p25': st_full['p25'],
            'full_asl_p75': st_full['p75'],
            'full_asl_p95': st_full['p95'],
            'full_asl_final10ns_mean': st_full['final_10ns_mean'],
            # Stem
            'stem_mean': st_stem['mean'],
            'stem_std': st_stem['std'],
            'stem_median': st_stem['median'],
            'stem_p25': st_stem['p25'],
            'stem_p75': st_stem['p75'],
            'stem_p95': st_stem['p95'],
            'stem_final10ns_mean': st_stem['final_10ns_mean'],
            # Loop Global (stem aligned)
            'loop_global_mean': st_loop_glob['mean'],
            'loop_global_std': st_loop_glob['std'],
            'loop_global_median': st_loop_glob['median'],
            'loop_global_p25': st_loop_glob['p25'],
            'loop_global_p75': st_loop_glob['p75'],
            'loop_global_p95': st_loop_glob['p95'],
            'loop_global_final10ns_mean': st_loop_glob['final_10ns_mean'],
            # Loop Internal (loop aligned)
            'loop_internal_mean': st_loop_int['mean'],
            'loop_internal_std': st_loop_int['std'],
            'loop_internal_median': st_loop_int['median'],
            'loop_internal_p25': st_loop_int['p25'],
            'loop_internal_p75': st_loop_int['p75'],
            'loop_internal_p95': st_loop_int['p95'],
            'loop_internal_final10ns_mean': st_loop_int['final_10ns_mean']
        }
        replicate_rows.append(row)
        
        print(f"[{rep_tag}]")
        print(f"  Full ASL RMSD: {st_full['mean']:.3f} +/- {st_full['std']:.3f} A (median: {st_full['median']:.3f} A)")
        print(f"  Stem RMSD:     {st_stem['mean']:.3f} +/- {st_stem['std']:.3f} A (median: {st_stem['median']:.3f} A)")
        print(f"  Loop Global:   {st_loop_glob['mean']:.3f} +/- {st_loop_glob['std']:.3f} A (median: {st_loop_glob['median']:.3f} A)")
        print(f"  Loop Internal: {st_loop_int['mean']:.3f} +/- {st_loop_int['std']:.3f} A (median: {st_loop_int['median']:.3f} A)")

    # Save Tables
    df_reps = pd.DataFrame(replicate_rows)
    df_reps.to_csv(os.path.join(TABLES_DIR, 'rmsd_replicates.csv'), index=False)
    print(f"\nSaved replicate-level RMSD table: {os.path.join(TABLES_DIR, 'rmsd_replicates.csv')}")
    
    # Save Time Series
    df_all_ts = pd.concat(time_series_dfs, ignore_index=True)
    df_all_ts.to_csv(os.path.join(TABLES_DIR, 'rmsd_time_series.csv'), index=False)
    print(f"Saved full time series table: {os.path.join(TABLES_DIR, 'rmsd_time_series.csv')}")
    
    # Construct-level summary (mean +/- SD of replicate means)
    summary_rows = []
    for c in ['wt', 'cua']:
        sub = df_reps[df_reps['construct'] == c]
        summary_rows.append({
            'construct': c,
            'n_replicates': len(sub),
            'full_asl_mean_of_means': float(sub['full_asl_mean'].mean()),
            'full_asl_std_of_means': float(sub['full_asl_mean'].std()),
            'stem_mean_of_means': float(sub['stem_mean'].mean()),
            'stem_std_of_means': float(sub['stem_mean'].std()),
            'loop_global_mean_of_means': float(sub['loop_global_mean'].mean()),
            'loop_global_std_of_means': float(sub['loop_global_mean'].std()),
            'loop_internal_mean_of_means': float(sub['loop_internal_mean'].mean()),
            'loop_internal_std_of_means': float(sub['loop_internal_mean'].std()),
            'loop_global_final10ns_mean_of_means': float(sub['loop_global_final10ns_mean'].mean()),
            'loop_global_final10ns_std_of_means': float(sub['loop_global_final10ns_mean'].std())
        })
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(os.path.join(TABLES_DIR, 'rmsd_summary.csv'), index=False)
    print(f"Saved construct-level RMSD summary: {os.path.join(TABLES_DIR, 'rmsd_summary.csv')}")
    print("=" * 80 + "\n")

if __name__ == '__main__':
    run_rmsd_analysis()
