"""
Phase 6: Anticodon Loop 3'-Stacking Geometry Analysis.
Analyzes base stacking across adjacent nucleotides in the 3' loop region:
  - Step 34-35 (anticodon base 1 to base 2)
  - Step 35-36 (anticodon base 2 to base 3)
  - Step 36-37 (anticodon base 3 to purine 37)
  - Step 37-38 (purine 37 to nucleotide 38)

Geometric criteria for nucleobase stacking:
  1. Base plane centroid distance R_centroid <= 4.5 Angstroms (0.45 nm)
  2. Inter-plane angle theta <= 30.0 degrees (cos theta >= 0.866)
  Ring atoms used for centroid and plane determination (via SVD best-fit plane):
    - Pyrimidines (C, U): N1, C2, N3, C4, C5, C6 (6-ring)
    - Purines (A, G): N1, C2, N3, C4, C5, C6, N7, C8, N9 (9-ring)

Outputs:
  - analysis/tables/stacking_replicates.csv
  - analysis/tables/stacking_summary.csv
  - analysis/tables/stacking_time_series.csv
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

# Steps along 3'-anticodon loop
# (step_name, posA, posB, res_idxA, res_idxB)
STACKING_STEPS = [
    ('34-35', 34, 35, 7, 8),
    ('35-36', 35, 36, 8, 9),
    ('36-37', 36, 37, 9, 10),
    ('37-38', 37, 38, 10, 11),
]

PYRIMIDINE_RINGS = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']
PURINE_RINGS = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9']

DIST_CUTOFF_STANDARD = 4.5   # Angstroms
ANGLE_CUTOFF_STANDARD = 30.0 # degrees

DIST_CUTOFF_PERMISSIVE = 5.0 # Angstroms
ANGLE_CUTOFF_PERMISSIVE = 35.0 # degrees

def get_base_ring_atoms(res):
    if res.name in ['A', 'G']:
        target_names = PURINE_RINGS
    else:
        target_names = PYRIMIDINE_RINGS
    return [a.index for a in res.atoms if a.name in target_names]

def compute_centroid_and_normal(traj_xyz, atom_indices):
    """
    traj_xyz: (n_frames, n_atoms, 3) in Angstroms
    atom_indices: list of ints
    Returns:
      centroids: (n_frames, 3)
      normals: (n_frames, 3) unit vectors
    """
    coords = traj_xyz[:, atom_indices, :] # (n_frames, k, 3)
    centroids = np.mean(coords, axis=1)   # (n_frames, 3)
    
    centered = coords - centroids[:, np.newaxis, :] # (n_frames, k, 3)
    
    # Compute SVD for best-fit plane normal (eigenvector with smallest eigenvalue)
    # Using cov matrix M = centered.T @ centered for each frame
    n_frames = coords.shape[0]
    normals = np.zeros((n_frames, 3))
    
    for f in range(n_frames):
        # centered[f] is (k, 3)
        _, _, vh = np.linalg.svd(centered[f], full_matrices=False)
        # vh[2] is normal to the best-fit plane
        normals[f] = vh[2]
        
    # Ensure unit length
    normals = normals / (np.linalg.norm(normals, axis=1, keepdims=True) + 1e-9)
    return centroids, normals

def run_stacking_analysis():
    print("=" * 80)
    print("PHASE 6: 3'-ANTICODON BASE STACKING GEOMETRY ANALYSIS")
    print("=" * 80)
    
    rep_records = []
    ts_list = []
    
    for construct, rep_num in REPLICATES:
        rep_tag = f"{construct}_rep{rep_num}"
        dcd_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_imaged.dcd')
        pdb_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_ref.pdb')
        
        traj = md.load(dcd_path, top=pdb_path)
        n_frames = traj.n_frames
        times_ns = np.linspace(0.010, 20.000, n_frames)
        xyz_ang = traj.xyz * 10.0 # Angstroms
        
        print(f"[{rep_tag}] Analyzing stacking across {n_frames} frames...")
        
        for step_name, posA, posB, rA, rB in STACKING_STEPS:
            resA = traj.topology.residue(rA)
            resB = traj.topology.residue(rB)
            
            atomsA = get_base_ring_atoms(resA)
            atomsB = get_base_ring_atoms(resB)
            
            centA, normA = compute_centroid_and_normal(xyz_ang, atomsA)
            centB, normB = compute_centroid_and_normal(xyz_ang, atomsB)
            
            # Centroid distance
            r_cent = np.linalg.norm(centA - centB, axis=-1)
            
            # Angle between plane normals (acute angle between non-oriented planes)
            cos_theta = np.abs(np.sum(normA * normB, axis=-1))
            cos_theta = np.clip(cos_theta, 0.0, 1.0)
            theta_deg = np.degrees(np.arccos(cos_theta))
            
            # Stacking boolean criteria
            is_stacked = (r_cent <= DIST_CUTOFF_STANDARD) & (theta_deg <= ANGLE_CUTOFF_STANDARD)
            is_stacked_perm = (r_cent <= DIST_CUTOFF_PERMISSIVE) & (theta_deg <= ANGLE_CUTOFF_PERMISSIVE)
            
            final_slice = slice(1000, 2000)
            
            rep_records.append({
                'construct': construct,
                'replicate': rep_num,
                'replicate_tag': rep_tag,
                'step': step_name,
                'pair_label': f"{resA.name}{posA}-{resB.name}{posB}",
                'posA': posA,
                'posB': posB,
                'resA': resA.name,
                'resB': resB.name,
                'stacking_occupancy': float(np.mean(is_stacked)),
                'stacking_permissive_occupancy': float(np.mean(is_stacked_perm)),
                'stacking_final10ns_occupancy': float(np.mean(is_stacked[final_slice])),
                'dist_mean_A': float(np.mean(r_cent)),
                'dist_std_A': float(np.std(r_cent)),
                'dist_median_A': float(np.median(r_cent)),
                'angle_mean_deg': float(np.mean(theta_deg)),
                'angle_std_deg': float(np.std(theta_deg)),
                'angle_median_deg': float(np.median(theta_deg))
            })
            
            df_ts = pd.DataFrame({
                'frame': np.arange(n_frames),
                'time_ns': np.round(times_ns, 4),
                'construct': construct,
                'replicate': rep_num,
                'replicate_tag': rep_tag,
                'step': step_name,
                'pair_label': f"{resA.name}{posA}-{resB.name}{posB}",
                'centroid_dist_A': np.round(r_cent, 4),
                'tilt_angle_deg': np.round(theta_deg, 4),
                'is_stacked': is_stacked.astype(int),
                'is_stacked_permissive': is_stacked_perm.astype(int)
            })
            ts_list.append(df_ts)
            
            print(f"  Step {resA.name}{posA}-{resB.name}{posB}: "
                  f"Stacking occ = {np.mean(is_stacked)*100:.2f}% | "
                  f"Mean R = {np.mean(r_cent):.2f} A | Mean theta = {np.mean(theta_deg):.1f} deg")

    df_rep = pd.DataFrame(rep_records)
    rep_csv = os.path.join(TABLES_DIR, 'stacking_replicates.csv')
    df_rep.to_csv(rep_csv, index=False)
    print(f"Saved replicate-level stacking table: {rep_csv}")
    
    df_ts_all = pd.concat(ts_list, ignore_index=True)
    ts_csv = os.path.join(TABLES_DIR, 'stacking_time_series.csv')
    df_ts_all.to_csv(ts_csv, index=False)
    print(f"Saved stacking time series table: {ts_csv}")
    
    # Construct-level summary across N=3 replicates
    summary_records = []
    for construct in ['wt', 'cua']:
        df_c = df_rep[df_rep['construct'] == construct]
        for step_name, posA, posB, _, _ in STACKING_STEPS:
            df_step = df_c[df_c['step'] == step_name]
            occ_vals = df_step['stacking_occupancy'].values
            occ_perm = df_step['stacking_permissive_occupancy'].values
            occ_final = df_step['stacking_final10ns_occupancy'].values
            dist_vals = df_step['dist_mean_A'].values
            ang_vals = df_step['angle_mean_deg'].values
            pair_lbl = df_step['pair_label'].iloc[0]
            
            summary_records.append({
                'construct': construct,
                'step': step_name,
                'pair_label': pair_lbl,
                'n_replicates': len(occ_vals),
                'stacking_occupancy_mean': float(np.mean(occ_vals)),
                'stacking_occupancy_std': float(np.std(occ_vals, ddof=1)),
                'stacking_occupancy_min': float(np.min(occ_vals)),
                'stacking_occupancy_max': float(np.max(occ_vals)),
                'stacking_permissive_mean': float(np.mean(occ_perm)),
                'stacking_permissive_std': float(np.std(occ_perm, ddof=1)),
                'stacking_final10ns_mean': float(np.mean(occ_final)),
                'stacking_final10ns_std': float(np.std(occ_final, ddof=1)),
                'centroid_dist_mean': float(np.mean(dist_vals)),
                'centroid_dist_std': float(np.std(dist_vals, ddof=1)),
                'tilt_angle_mean': float(np.mean(ang_vals)),
                'tilt_angle_std': float(np.std(ang_vals, ddof=1)),
            })
            
    df_summary = pd.DataFrame(summary_records)
    sum_csv = os.path.join(TABLES_DIR, 'stacking_summary.csv')
    df_summary.to_csv(sum_csv, index=False)
    print(f"Saved construct-level stacking summary: {sum_csv}")
    print("=" * 80)

if __name__ == '__main__':
    run_stacking_analysis()
