"""
Phase 6: Per-Residue Root-Mean-Square Fluctuation (RMSF) Analysis.
Calculates per-residue RMSF for all 17 ASL nucleotides across all 6 replicates:
  - WT (replicates 1, 2, 3)
  - CUA (replicates 1, 2, 3)

Alignment reference: Stem backbone heavy atoms (P, O5', C5', C4', C3', O3' of residues 27-31 and 39-43)
of the Phase 4C final equilibrated structure. This measures local residue fluctuations as well
as loop mobility relative to the stable stem scaffold.

Outputs:
  - analysis/tables/rmsf_replicates.csv
  - analysis/tables/rmsf_summary.csv
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

STEM_RES_INDICES = [0, 1, 2, 3, 4, 12, 13, 14, 15, 16] # tRNA 27-31, 39-43
LOOP_RES_INDICES = [5, 6, 7, 8, 9, 10, 11]              # tRNA 32-38
BACKBONE_ATOMS = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]

TRNA_POSITIONS = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43]

def run_rmsf_analysis():
    print("=" * 80)
    print("PHASE 6: PER-RESIDUE RMSF ANALYSIS")
    print("=" * 80)
    
    rep_records = []
    
    for construct, rep_num in REPLICATES:
        rep_tag = f"{construct}_rep{rep_num}"
        dcd_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_imaged.dcd')
        pdb_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_ref.pdb')
        
        # Load reference
        p4c_pdb = os.path.join(BASE_DIR, 'simulations', construct, 'equilibration', 'final_equilibrated.pdb')
        p4c_full = md.load(p4c_pdb)
        solute_idx = p4c_full.topology.select('not water and not (resname K or resname CL)')
        ref = p4c_full.atom_slice(solute_idx)
        
        traj = md.load(dcd_path, top=pdb_path)
        
        # Identify stem backbone atoms for superposition
        stem_bb = [
            a.index for a in traj.topology.atoms
            if a.name in BACKBONE_ATOMS and a.residue.index in STEM_RES_INDICES
        ]
        
        # Align trajectory to reference stem backbone
        traj.superpose(ref, atom_indices=stem_bb)
        
        # Compute atomic RMSF (in nm, convert to Angstrom)
        # md.rmsf computes fluctuation around average structure of aligned trajectory
        atomic_rmsf = md.rmsf(traj, traj, atom_indices=None) * 10.0 # Angstroms
        
        for r_idx in range(traj.topology.n_residues):
            res = traj.topology.residue(r_idx)
            trna_pos = TRNA_POSITIONS[r_idx]
            region = "stem" if r_idx in STEM_RES_INDICES else "loop"
            
            # Backbone atoms of this residue
            res_bb_atoms = [a.index for a in res.atoms if a.name in BACKBONE_ATOMS]
            # Heavy atoms of this residue (non-hydrogen)
            res_heavy_atoms = [a.index for a in res.atoms if a.element.symbol != 'H']
            
            # Per-residue RMSF is defined as root-mean-square of atomic fluctuations in residue
            res_bb_rmsf = np.sqrt(np.mean(atomic_rmsf[res_bb_atoms]**2)) if res_bb_atoms else 0.0
            res_heavy_rmsf = np.sqrt(np.mean(atomic_rmsf[res_heavy_atoms]**2)) if res_heavy_atoms else 0.0
            
            rep_records.append({
                'construct': construct,
                'replicate': rep_num,
                'replicate_tag': rep_tag,
                'residue_idx_0based': r_idx,
                'residue_idx_1based': r_idx + 1,
                'trna_position': trna_pos,
                'resname': res.name,
                'region': region,
                'backbone_rmsf_angstrom': float(res_bb_rmsf),
                'heavy_atom_rmsf_angstrom': float(res_heavy_rmsf)
            })
            
        print(f"[{rep_tag}] RMSF calculated across {traj.topology.n_residues} residues.")

    df_rep = pd.DataFrame(rep_records)
    rep_csv_path = os.path.join(TABLES_DIR, 'rmsf_replicates.csv')
    df_rep.to_csv(rep_csv_path, index=False)
    print(f"Saved replicate-level RMSF table: {rep_csv_path}")
    
    # Construct-level summary (mean ± SD across the N=3 independent replicates)
    summary_rows = []
    for construct in ['wt', 'cua']:
        df_c = df_rep[df_rep['construct'] == construct]
        for trna_pos in TRNA_POSITIONS:
            df_pos = df_c[df_c['trna_position'] == trna_pos]
            r_idx_0 = df_pos['residue_idx_0based'].iloc[0]
            resname = df_pos['resname'].iloc[0]
            region = df_pos['region'].iloc[0]
            
            bb_vals = df_pos['backbone_rmsf_angstrom'].values
            heavy_vals = df_pos['heavy_atom_rmsf_angstrom'].values
            
            summary_rows.append({
                'construct': construct,
                'trna_position': trna_pos,
                'residue_idx_1based': r_idx_0 + 1,
                'resname': resname,
                'region': region,
                'n_replicates': len(bb_vals),
                'backbone_rmsf_mean': float(np.mean(bb_vals)),
                'backbone_rmsf_std': float(np.std(bb_vals, ddof=1)),
                'backbone_rmsf_min': float(np.min(bb_vals)),
                'backbone_rmsf_max': float(np.max(bb_vals)),
                'heavy_rmsf_mean': float(np.mean(heavy_vals)),
                'heavy_rmsf_std': float(np.std(heavy_vals, ddof=1)),
                'heavy_rmsf_min': float(np.min(heavy_vals)),
                'heavy_rmsf_max': float(np.max(heavy_vals))
            })
            
    df_summary = pd.DataFrame(summary_rows)
    summary_csv_path = os.path.join(TABLES_DIR, 'rmsf_summary.csv')
    df_summary.to_csv(summary_csv_path, index=False)
    print(f"Saved construct-level RMSF summary: {summary_csv_path}")
    print("=" * 80)

if __name__ == '__main__':
    run_rmsf_analysis()
