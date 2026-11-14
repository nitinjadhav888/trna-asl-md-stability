"""
Phase 6: Stem Base-Pair Hydrogen-Bond Occupancy Analysis.
Analyzes the 5 canonical Watson-Crick base pairs in the ASL stem:
  Pair 1: C27 - G43 (residue 0 - residue 16)
  Pair 2: C28 - G42 (residue 1 - residue 15)
  Pair 3: C29 - G41 (residue 2 - residue 14)
  Pair 4: G30 - C40 (residue 3 - residue 13)
  Pair 5: C31 - G39 (residue 4 - residue 12)

Standard geometric hydrogen-bond criteria:
  - Donor heavy atom to Acceptor heavy atom distance <= 3.5 Angstroms (0.35 nm)
  - Hydrogen atom to Acceptor heavy atom distance <= 2.5 Angstroms (0.25 nm)
  - Donor - H ... Acceptor angle >= 120.0 degrees

Watson-Crick G-C hydrogen bonds:
  HB1: G(N1)-H1 ... C(N3)
  HB2: G(N2)-H21/H22 ... C(O2)
  HB3: C(N4)-H41/H42 ... G(O6)

Base pair intact criterion: >= 2 of 3 canonical hydrogen bonds formed simultaneously.

Outputs:
  - analysis/tables/stem_hbond_replicates.csv
  - analysis/tables/stem_hbond_summary.csv
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

# Canonical 5 base pairs: (name, res_idx_5p, res_idx_3p, type)
# where 5p is 5' strand, 3p is 3' strand
STEM_PAIRS = [
    ('C27:G43', 0, 16, 'C:G'),
    ('C28:G42', 1, 15, 'C:G'),
    ('C29:G41', 2, 14, 'C:G'),
    ('G30:C40', 3, 13, 'G:C'),
    ('C31:G39', 4, 12, 'C:G'),
]

DIST_HEAVY_CUTOFF = 0.35 # nm (3.5 A)
DIST_H_CUTOFF = 0.25     # nm (2.5 A)
ANGLE_CUTOFF = 120.0     # degrees

def compute_angle(xyz, idx_d, idx_h, idx_a):
    """
    Computes angle D-H...A in degrees for all frames.
    xyz: (n_frames, n_atoms, 3)
    """
    v_hd = xyz[:, idx_d, :] - xyz[:, idx_h, :] # H -> D
    v_ha = xyz[:, idx_a, :] - xyz[:, idx_h, :] # H -> A
    
    norm_hd = np.linalg.norm(v_hd, axis=-1)
    norm_ha = np.linalg.norm(v_ha, axis=-1)
    
    dot = np.sum(v_hd * v_ha, axis=-1)
    cosine = np.clip(dot / (norm_hd * norm_ha + 1e-9), -1.0, 1.0)
    angle_rad = np.arccos(cosine)
    return np.degrees(angle_rad)

def eval_hbond(traj, d_name, h_names, a_name, res_d_idx, res_a_idx):
    """
    Evaluates boolean array (n_frames,) indicating whether D-H...A HB is formed.
    If multiple H atoms (e.g. amino group), HB formed if ANY satisfies criterion.
    """
    res_d = traj.topology.residue(res_d_idx)
    res_a = traj.topology.residue(res_a_idx)
    
    d_atom = [a.index for a in res_d.atoms if a.name == d_name][0]
    a_atom = [a.index for a in res_a.atoms if a.name == a_name][0]
    h_atoms = [a.index for a in res_d.atoms if a.name in h_names]
    
    # Heavy atom distance (n_frames,)
    d_da = np.linalg.norm(traj.xyz[:, d_atom, :] - traj.xyz[:, a_atom, :], axis=-1)
    
    hb_formed = np.zeros(traj.n_frames, dtype=bool)
    
    for h_atom in h_atoms:
        d_ha = np.linalg.norm(traj.xyz[:, h_atom, :] - traj.xyz[:, a_atom, :], axis=-1)
        angles = compute_angle(traj.xyz, d_atom, h_atom, a_atom)
        h_satisfies = (d_da <= DIST_HEAVY_CUTOFF) & (d_ha <= DIST_H_CUTOFF) & (angles >= ANGLE_CUTOFF)
        hb_formed = hb_formed | h_satisfies
        
    return hb_formed, d_da * 10.0 # return boolean array and D-A distance in Angstroms

def run_hbond_analysis():
    print("=" * 80)
    print("PHASE 6: STEM BASE-PAIR HYDROGEN-BOND OCCUPANCY ANALYSIS")
    print("=" * 80)
    
    replicate_records = []
    
    for construct, rep_num in REPLICATES:
        rep_tag = f"{construct}_rep{rep_num}"
        dcd_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_imaged.dcd')
        pdb_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_ref.pdb')
        
        traj = md.load(dcd_path, top=pdb_path)
        n_frames = traj.n_frames
        
        print(f"[{rep_tag}] Evaluating 5 stem base pairs across {n_frames} frames...")
        
        for pair_name, r5, r3, bp_type in STEM_PAIRS:
            # Determine which residue is C and which is G
            if bp_type == 'C:G':
                c_idx, g_idx = r5, r3
            else:
                g_idx, c_idx = r5, r3
                
            # HB1: G(N1)-H1 ... C(N3)
            hb1_arr, hb1_dist = eval_hbond(traj, 'N1', ['H1'], 'N3', g_idx, c_idx)
            # HB2: G(N2)-H21/22 ... C(O2)
            hb2_arr, hb2_dist = eval_hbond(traj, 'N2', ['H21', 'H22'], 'O2', g_idx, c_idx)
            # HB3: C(N4)-H41/42 ... G(O6)
            hb3_arr, hb3_dist = eval_hbond(traj, 'N4', ['H41', 'H42'], 'O6', c_idx, g_idx)
            
            # Intact if >= 2 H-bonds present
            intact_arr = (hb1_arr.astype(int) + hb2_arr.astype(int) + hb3_arr.astype(int)) >= 2
            
            # Final 10 ns slice
            final_slice = slice(1000, 2000)
            
            replicate_records.append({
                'construct': construct,
                'replicate': rep_num,
                'replicate_tag': rep_tag,
                'base_pair': pair_name,
                'pair_type': 'Watson-Crick G-C',
                'res_idx_5p': r5,
                'res_idx_3p': r3,
                'hb1_occupancy': float(np.mean(hb1_arr)),
                'hb2_occupancy': float(np.mean(hb2_arr)),
                'hb3_occupancy': float(np.mean(hb3_arr)),
                'bp_intact_occupancy': float(np.mean(intact_arr)),
                'bp_intact_final10ns': float(np.mean(intact_arr[final_slice])),
                'hb1_mean_dist_A': float(np.mean(hb1_dist)),
                'hb2_mean_dist_A': float(np.mean(hb2_dist)),
                'hb3_mean_dist_A': float(np.mean(hb3_dist)),
            })
            
    df_rep = pd.DataFrame(replicate_records)
    rep_csv = os.path.join(TABLES_DIR, 'stem_hbond_replicates.csv')
    df_rep.to_csv(rep_csv, index=False)
    print(f"Saved replicate-level stem H-bond table: {rep_csv}")
    
    # Construct-level summary across N=3 replicates
    summary_records = []
    for construct in ['wt', 'cua']:
        df_c = df_rep[df_rep['construct'] == construct]
        for pair_name, _, _, _ in STEM_PAIRS:
            df_pair = df_c[df_c['base_pair'] == pair_name]
            
            intact_vals = df_pair['bp_intact_occupancy'].values
            final_vals = df_pair['bp_intact_final10ns'].values
            hb1_vals = df_pair['hb1_occupancy'].values
            hb2_vals = df_pair['hb2_occupancy'].values
            hb3_vals = df_pair['hb3_occupancy'].values
            
            summary_records.append({
                'construct': construct,
                'base_pair': pair_name,
                'n_replicates': len(intact_vals),
                'bp_intact_mean': float(np.mean(intact_vals)),
                'bp_intact_std': float(np.std(intact_vals, ddof=1)),
                'bp_intact_min': float(np.min(intact_vals)),
                'bp_intact_max': float(np.max(intact_vals)),
                'bp_intact_final10ns_mean': float(np.mean(final_vals)),
                'bp_intact_final10ns_std': float(np.std(final_vals, ddof=1)),
                'hb1_mean_occupancy': float(np.mean(hb1_vals)),
                'hb2_mean_occupancy': float(np.mean(hb2_vals)),
                'hb3_mean_occupancy': float(np.mean(hb3_vals)),
            })
            
    df_summary = pd.DataFrame(summary_records)
    sum_csv = os.path.join(TABLES_DIR, 'stem_hbond_summary.csv')
    df_summary.to_csv(sum_csv, index=False)
    print(f"Saved construct-level stem H-bond summary: {sum_csv}")
    print("=" * 80)

if __name__ == '__main__':
    run_hbond_analysis()
