"""
Phase 6: Additional Structural Metric — Sugar Pucker at Nucleotide 34.
Calculates the ribose pseudorotation phase angle (P), maximum pucker amplitude (tau_m),
and conformational states (North/C3'-endo vs South/C2'-endo) for nucleotide 34
(G34 in WT, C34 in CUA) across all 6 production replicates.

Theoretical formulation (Altona & Sundaralingam, JACS 1972):
  Endocyclic furanose ring torsions:
    nu0: C4' - O4' - C1' - C2'
    nu1: O4' - C1' - C2' - C3'
    nu2: C1' - C2' - C3' - C4'
    nu3: C2' - C3' - C4' - O4'
    nu4: C3' - C4' - O4' - C1'
  Pseudorotation phase P:
    tan(P) = [ (nu4 + nu1) - (nu3 + nu0) ] / [ 2 * nu2 * (sin 36 deg + sin 72 deg) ]
    P in [0, 360) degrees
  Maximum amplitude tau_m:
    tau_m = nu2 / cos(P)

Conformational classification:
  - North (C3'-endo, A-form canonical RNA): 0 <= P < 36 deg or 324 <= P < 360 deg
  - South (C2'-endo, B-form / extended loop): 144 <= P < 216 deg
  - East (O4'-endo): 36 <= P < 144 deg
  - West (C1'-exo / C4'-endo): 216 <= P < 324 deg

Outputs:
  - analysis/tables/sugar_pucker_replicates.csv
  - analysis/tables/sugar_pucker_summary.csv
  - analysis/tables/sugar_pucker_time_series.csv
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

POS34_RES_IDX = 7 # Residue 8 in 1-based, tRNA pos 34 (G in WT, C in CUA)

SIN_FACTOR = np.sin(np.radians(36.0)) + np.sin(np.radians(72.0)) # ~ 1.5388417685876268

def compute_pseudorotation(nu0, nu1, nu2, nu3, nu4):
    """
    nu0..nu4 are arrays of shape (n_frames,) in degrees.
    Returns P in degrees [0, 360) and tau_m in degrees.
    """
    numerator = (nu4 + nu1) - (nu3 + nu0)
    denominator = 2.0 * nu2 * SIN_FACTOR
    
    P_rad = np.arctan2(numerator, denominator)
    P_deg = np.degrees(P_rad)
    P_deg = np.where(P_deg < 0, P_deg + 360.0, P_deg)
    
    # Amplitude tau_m
    cos_P = np.cos(np.radians(P_deg))
    tau_m = np.where(np.abs(cos_P) > 1e-4, nu2 / cos_P, np.abs(nu2))
    tau_m = np.abs(tau_m)
    
    return P_deg, tau_m

def classify_pucker(P_deg):
    """
    Classifies pseudorotation phase into conformational state.
    """
    state = np.full(P_deg.shape, 'Other', dtype=object)
    
    # North (C3'-endo): 0 <= P < 36 or 324 <= P < 360
    is_north = (P_deg < 36.0) | (P_deg >= 324.0)
    # South (C2'-endo): 144 <= P < 216
    is_south = (P_deg >= 144.0) & (P_deg < 216.0)
    # East (O4'-endo): 36 <= P < 144
    is_east = (P_deg >= 36.0) & (P_deg < 144.0)
    # West (C1'-exo / C4'-endo): 216 <= P < 324
    is_west = (P_deg >= 216.0) & (P_deg < 324.0)
    
    state[is_north] = 'North (C3\'-endo)'
    state[is_south] = 'South (C2\'-endo)'
    state[is_east] = 'East (O4\'-endo)'
    state[is_west] = 'West (C1\'-exo/C4\'-endo)'
    
    return state, is_north, is_south

def run_pucker_analysis():
    print("=" * 80)
    print("PHASE 6: RIBOSE SUGAR PUCKER ANALYSIS AT NUCLEOTIDE 34")
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
        
        res = traj.topology.residue(POS34_RES_IDX)
        
        # Atom map for ribose ring
        c1p = [a.index for a in res.atoms if a.name == "C1'"][0]
        c2p = [a.index for a in res.atoms if a.name == "C2'"][0]
        c3p = [a.index for a in res.atoms if a.name == "C3'"][0]
        c4p = [a.index for a in res.atoms if a.name == "C4'"][0]
        o4p = [a.index for a in res.atoms if a.name == "O4'"][0]
        
        # Compute the 5 dihedrals (radians to degrees)
        dihed_atoms = np.array([
            [c4p, o4p, c1p, c2p], # nu0
            [o4p, c1p, c2p, c3p], # nu1
            [c1p, c2p, c3p, c4p], # nu2
            [c2p, c3p, c4p, o4p], # nu3
            [c3p, c4p, o4p, c1p], # nu4
        ])
        
        diheds_rad = md.compute_dihedrals(traj, dihed_atoms) # (n_frames, 5)
        diheds_deg = np.degrees(diheds_rad)
        
        nu0 = diheds_deg[:, 0]
        nu1 = diheds_deg[:, 1]
        nu2 = diheds_deg[:, 2]
        nu3 = diheds_deg[:, 3]
        nu4 = diheds_deg[:, 4]
        
        P_deg, tau_m = compute_pseudorotation(nu0, nu1, nu2, nu3, nu4)
        states, is_north, is_south = classify_pucker(P_deg)
        
        final_slice = slice(1000, 2000)
        
        north_occ = float(np.mean(is_north))
        south_occ = float(np.mean(is_south))
        east_occ = float(np.mean(states == 'East (O4\'-endo)'))
        west_occ = float(np.mean(states == 'West (C1\'-exo/C4\'-endo)'))
        
        rep_records.append({
            'construct': construct,
            'replicate': rep_num,
            'replicate_tag': rep_tag,
            'target_res': f"{res.name}34",
            'frames': n_frames,
            'north_c3_endo_occupancy': north_occ,
            'south_c2_endo_occupancy': south_occ,
            'east_o4_endo_occupancy': east_occ,
            'west_c1_exo_occupancy': west_occ,
            'north_final10ns_occupancy': float(np.mean(is_north[final_slice])),
            'south_final10ns_occupancy': float(np.mean(is_south[final_slice])),
            'mean_phase_deg': float(np.mean(P_deg)),
            'std_phase_deg': float(np.std(P_deg)),
            'median_phase_deg': float(np.median(P_deg)),
            'mean_amplitude_deg': float(np.mean(tau_m)),
            'std_amplitude_deg': float(np.std(tau_m))
        })
        
        df_ts = pd.DataFrame({
            'frame': np.arange(n_frames),
            'time_ns': np.round(times_ns, 4),
            'construct': construct,
            'replicate': rep_num,
            'replicate_tag': rep_tag,
            'target_res': f"{res.name}34",
            'nu0_deg': np.round(nu0, 2),
            'nu1_deg': np.round(nu1, 2),
            'nu2_deg': np.round(nu2, 2),
            'nu3_deg': np.round(nu3, 2),
            'nu4_deg': np.round(nu4, 2),
            'phase_P_deg': np.round(P_deg, 2),
            'amplitude_tau_m_deg': np.round(tau_m, 2),
            'conformational_state': states,
            'is_north': is_north.astype(int),
            'is_south': is_south.astype(int)
        })
        ts_list.append(df_ts)
        
        print(f"[{rep_tag}] {res.name}34 Pucker: North (C3'-endo) = {north_occ*100:.2f}% | "
              f"South (C2'-endo) = {south_occ*100:.2f}% | East = {east_occ*100:.2f}% | West = {west_occ*100:.2f}% | "
              f"Median P = {np.median(P_deg):.1f} deg | Mean tau_m = {np.mean(tau_m):.1f} deg")

    df_rep = pd.DataFrame(rep_records)
    rep_csv = os.path.join(TABLES_DIR, 'sugar_pucker_replicates.csv')
    df_rep.to_csv(rep_csv, index=False)
    print(f"Saved replicate-level sugar pucker table: {rep_csv}")
    
    df_ts_all = pd.concat(ts_list, ignore_index=True)
    ts_csv = os.path.join(TABLES_DIR, 'sugar_pucker_time_series.csv')
    df_ts_all.to_csv(ts_csv, index=False)
    print(f"Saved sugar pucker time series table: {ts_csv}")
    
    # Construct-level summary across N=3 replicates
    summary_records = []
    for construct in ['wt', 'cua']:
        df_c = df_rep[df_rep['construct'] == construct]
        n_vals = df_c['north_c3_endo_occupancy'].values
        s_vals = df_c['south_c2_endo_occupancy'].values
        e_vals = df_c['east_o4_endo_occupancy'].values
        w_vals = df_c['west_c1_exo_occupancy'].values
        n_final_vals = df_c['north_final10ns_occupancy'].values
        s_final_vals = df_c['south_final10ns_occupancy'].values
        p_med_vals = df_c['median_phase_deg'].values
        amp_vals = df_c['mean_amplitude_deg'].values
        target_res = df_c['target_res'].iloc[0]
        
        summary_records.append({
            'construct': construct,
            'target_res': target_res,
            'n_replicates': len(n_vals),
            'north_c3_endo_mean': float(np.mean(n_vals)),
            'north_c3_endo_std': float(np.std(n_vals, ddof=1)),
            'north_c3_endo_min': float(np.min(n_vals)),
            'north_c3_endo_max': float(np.max(n_vals)),
            'south_c2_endo_mean': float(np.mean(s_vals)),
            'south_c2_endo_std': float(np.std(s_vals, ddof=1)),
            'south_c2_endo_min': float(np.min(s_vals)),
            'south_c2_endo_max': float(np.max(s_vals)),
            'east_o4_endo_mean': float(np.mean(e_vals)),
            'east_o4_endo_std': float(np.std(e_vals, ddof=1)),
            'west_c1_exo_mean': float(np.mean(w_vals)),
            'west_c1_exo_std': float(np.std(w_vals, ddof=1)),
            'north_final10ns_mean': float(np.mean(n_final_vals)),
            'north_final10ns_std': float(np.std(n_final_vals, ddof=1)),
            'south_final10ns_mean': float(np.mean(s_final_vals)),
            'south_final10ns_std': float(np.std(s_final_vals, ddof=1)),
            'median_phase_deg_mean': float(np.mean(p_med_vals)),
            'amplitude_tau_m_mean': float(np.mean(amp_vals)),
            'amplitude_tau_m_std': float(np.std(amp_vals, ddof=1)),
        })
        
    df_summary = pd.DataFrame(summary_records)
    sum_csv = os.path.join(TABLES_DIR, 'sugar_pucker_summary.csv')
    df_summary.to_csv(sum_csv, index=False)
    print(f"Saved construct-level sugar pucker summary: {sum_csv}")
    print("=" * 80)

if __name__ == '__main__':
    run_pucker_analysis()
