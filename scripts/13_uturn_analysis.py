"""
Phase 6: Canonical tRNA Loop U-Turn Diagnostic Analysis.
Analyzes the U-turn diagnostic interaction:
  U33(N3) ... phosphate of position 36
  - In WT: position 36 is C36 (residue index 9)
  - In CUA: position 36 is A36 (residue index 9, suppressor mutation C36A)

Specific atoms evaluated:
  - Donor: U33(N3) (with proton H3)
  - Acceptors: non-bridging phosphate oxygens OP1 and OP2 of position 36
  - Phosphorus atom: P of position 36

Criteria:
  1. H-bond Criterion:
     - Heavy atom distance min(d(N3...OP1), d(N3...OP2)) <= 3.5 Angstroms
     - Proton distance min(d(H3...OP1), d(H3...OP2)) <= 2.5 Angstroms
     - Angle N3-H3...OP >= 120.0 degrees
  2. Contact Criterion:
     - Heavy atom distance min(d(N3...OP1), d(N3...OP2)) <= 3.5 Angstroms
     - Lenient contact cutoff: <= 4.0 Angstroms

Outputs:
  - analysis/tables/uturn_replicates.csv
  - analysis/tables/uturn_summary.csv
  - analysis/tables/uturn_time_series.csv
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

U33_RES_IDX = 6 # Residue 7 in 1-based, tRNA pos 33
POS36_RES_IDX = 9 # Residue 10 in 1-based, tRNA pos 36 (C in WT, A in CUA)

DIST_HEAVY_CUTOFF = 0.35 # nm (3.5 A)
DIST_H_CUTOFF = 0.25     # nm (2.5 A)
ANGLE_CUTOFF = 120.0     # degrees

def compute_angle(xyz, idx_d, idx_h, idx_a):
    v_hd = xyz[:, idx_d, :] - xyz[:, idx_h, :]
    v_ha = xyz[:, idx_a, :] - xyz[:, idx_h, :]
    norm_hd = np.linalg.norm(v_hd, axis=-1)
    norm_ha = np.linalg.norm(v_ha, axis=-1)
    dot = np.sum(v_hd * v_ha, axis=-1)
    cosine = np.clip(dot / (norm_hd * norm_ha + 1e-9), -1.0, 1.0)
    return np.degrees(np.arccos(cosine))

def run_uturn_analysis():
    print("=" * 80)
    print("PHASE 6: U-TURN DIAGNOSTIC ANALYSIS (U33 N3 ... P36 PHOSPHATE)")
    print("=" * 80)
    
    replicate_records = []
    ts_list = []
    
    for construct, rep_num in REPLICATES:
        rep_tag = f"{construct}_rep{rep_num}"
        dcd_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_imaged.dcd')
        pdb_path = os.path.join(BASE_DIR, 'simulations', rep_tag, 'trajectory', 'solute_ref.pdb')
        
        traj = md.load(dcd_path, top=pdb_path)
        n_frames = traj.n_frames
        times_ns = np.linspace(0.010, 20.000, n_frames)
        
        # Locate atoms
        u33_res = traj.topology.residue(U33_RES_IDX)
        p36_res = traj.topology.residue(POS36_RES_IDX)
        
        u33_n3 = [a.index for a in u33_res.atoms if a.name == 'N3'][0]
        u33_h3 = [a.index for a in u33_res.atoms if a.name == 'H3'][0]
        
        p36_p = [a.index for a in p36_res.atoms if a.name == 'P'][0]
        p36_op1 = [a.index for a in p36_res.atoms if a.name == 'OP1'][0]
        p36_op2 = [a.index for a in p36_res.atoms if a.name == 'OP2'][0]
        
        # Distances in nm
        d_n3_p = np.linalg.norm(traj.xyz[:, u33_n3, :] - traj.xyz[:, p36_p, :], axis=-1) * 10.0 # Angstroms
        d_n3_op1 = np.linalg.norm(traj.xyz[:, u33_n3, :] - traj.xyz[:, p36_op1, :], axis=-1) * 10.0
        d_n3_op2 = np.linalg.norm(traj.xyz[:, u33_n3, :] - traj.xyz[:, p36_op2, :], axis=-1) * 10.0
        d_min_op = np.minimum(d_n3_op1, d_n3_op2)
        
        d_h3_op1 = np.linalg.norm(traj.xyz[:, u33_h3, :] - traj.xyz[:, p36_op1, :], axis=-1) * 10.0
        d_h3_op2 = np.linalg.norm(traj.xyz[:, u33_h3, :] - traj.xyz[:, p36_op2, :], axis=-1) * 10.0
        
        ang_op1 = compute_angle(traj.xyz, u33_n3, u33_h3, p36_op1)
        ang_op2 = compute_angle(traj.xyz, u33_n3, u33_h3, p36_op2)
        
        # Criteria
        hb_op1 = (d_n3_op1 <= 3.5) & (d_h3_op1 <= 2.5) & (ang_op1 >= 120.0)
        hb_op2 = (d_n3_op2 <= 3.5) & (d_h3_op2 <= 2.5) & (ang_op2 >= 120.0)
        hb_intact = hb_op1 | hb_op2
        
        contact_35A = d_min_op <= 3.5
        contact_40A = d_min_op <= 4.0
        
        final_slice = slice(1000, 2000)
        
        # Summary row
        replicate_records.append({
            'construct': construct,
            'replicate': rep_num,
            'replicate_tag': rep_tag,
            'target_resname': p36_res.name,
            'hb_occupancy': float(np.mean(hb_intact)),
            'contact_35A_occupancy': float(np.mean(contact_35A)),
            'contact_40A_occupancy': float(np.mean(contact_40A)),
            'hb_final10ns_occupancy': float(np.mean(hb_intact[final_slice])),
            'contact_35A_final10ns': float(np.mean(contact_35A[final_slice])),
            'dist_min_mean_A': float(np.mean(d_min_op)),
            'dist_min_std_A': float(np.std(d_min_op)),
            'dist_min_median_A': float(np.median(d_min_op)),
            'dist_min_p25_A': float(np.percentile(d_min_op, 25)),
            'dist_min_p75_A': float(np.percentile(d_min_op, 75)),
            'dist_min_p95_A': float(np.percentile(d_min_op, 95)),
            'dist_p_mean_A': float(np.mean(d_n3_p)),
            'dist_p_median_A': float(np.median(d_n3_p)),
        })
        
        df_ts = pd.DataFrame({
            'frame': np.arange(n_frames),
            'time_ns': np.round(times_ns, 4),
            'construct': construct,
            'replicate': rep_num,
            'replicate_tag': rep_tag,
            'dist_n3_op_min_A': np.round(d_min_op, 4),
            'dist_n3_p_A': np.round(d_n3_p, 4),
            'hb_intact': hb_intact.astype(int),
            'contact_35A': contact_35A.astype(int),
            'contact_40A': contact_40A.astype(int),
        })
        ts_list.append(df_ts)
        
        print(f"[{rep_tag}] Target: {p36_res.name}36 | HB occ: {np.mean(hb_intact)*100:.2f}% | "
              f"Contact (<=3.5 A): {np.mean(contact_35A)*100:.2f}% | Mean dist: {np.mean(d_min_op):.2f} +/- {np.std(d_min_op):.2f} A")

    df_rep = pd.DataFrame(replicate_records)
    rep_csv = os.path.join(TABLES_DIR, 'uturn_replicates.csv')
    df_rep.to_csv(rep_csv, index=False)
    print(f"Saved replicate-level U-turn table: {rep_csv}")
    
    df_ts_all = pd.concat(ts_list, ignore_index=True)
    ts_csv = os.path.join(TABLES_DIR, 'uturn_time_series.csv')
    df_ts_all.to_csv(ts_csv, index=False)
    print(f"Saved U-turn time-series table: {ts_csv}")
    
    # Construct-level summary across N=3 replicates
    summary_rows = []
    for construct in ['wt', 'cua']:
        df_c = df_rep[df_rep['construct'] == construct]
        hb_vals = df_c['hb_occupancy'].values
        c35_vals = df_c['contact_35A_occupancy'].values
        c40_vals = df_c['contact_40A_occupancy'].values
        final_hb_vals = df_c['hb_final10ns_occupancy'].values
        dist_vals = df_c['dist_min_mean_A'].values
        dist_p_vals = df_c['dist_p_mean_A'].values
        
        summary_rows.append({
            'construct': construct,
            'target_res': 'C36' if construct == 'wt' else 'A36',
            'n_replicates': len(hb_vals),
            'hb_occupancy_mean': float(np.mean(hb_vals)),
            'hb_occupancy_std': float(np.std(hb_vals, ddof=1)),
            'hb_occupancy_min': float(np.min(hb_vals)),
            'hb_occupancy_max': float(np.max(hb_vals)),
            'contact_35A_mean': float(np.mean(c35_vals)),
            'contact_35A_std': float(np.std(c35_vals, ddof=1)),
            'contact_40A_mean': float(np.mean(c40_vals)),
            'contact_40A_std': float(np.std(c40_vals, ddof=1)),
            'hb_final10ns_mean': float(np.mean(final_hb_vals)),
            'hb_final10ns_std': float(np.std(final_hb_vals, ddof=1)),
            'dist_min_mean': float(np.mean(dist_vals)),
            'dist_min_std': float(np.std(dist_vals, ddof=1)),
            'dist_p_mean': float(np.mean(dist_p_vals)),
            'dist_p_std': float(np.std(dist_p_vals, ddof=1)),
        })
        
    df_summary = pd.DataFrame(summary_rows)
    sum_csv = os.path.join(TABLES_DIR, 'uturn_summary.csv')
    df_summary.to_csv(sum_csv, index=False)
    print(f"Saved construct-level U-turn summary: {sum_csv}")
    print("=" * 80)

if __name__ == '__main__':
    run_uturn_analysis()
