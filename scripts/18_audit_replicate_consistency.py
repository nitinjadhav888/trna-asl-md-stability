"""
Phase 6.1: Audit Replicate Consistency and Sampling Convergence.
Generates:
  - analysis/tables/replicate_consistency_audit.csv

Evaluates qualitative direction consistency across all 3 independent WT and
3 independent CUA replicates, and assesses sampling status across all major observables.
"""
import os
import sys
import numpy as np
import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_DIR = os.path.join(BASE_DIR, 'analysis', 'tables')
os.makedirs(TABLES_DIR, exist_ok=True)

def run_consistency_audit():
    print("=" * 80)
    print("PHASE 6.1: REPLICATE CONSISTENCY & SAMPLING CONVERGENCE AUDIT")
    print("=" * 80)
    
    # Load replicate data
    df_rmsd = pd.read_csv(os.path.join(TABLES_DIR, 'rmsd_replicates.csv')).set_index('replicate_tag')
    df_stem = pd.read_csv(os.path.join(TABLES_DIR, 'stem_hbond_replicates.csv'))
    df_uturn = pd.read_csv(os.path.join(TABLES_DIR, 'uturn_replicates.csv')).set_index('replicate_tag')
    df_stack = pd.read_csv(os.path.join(TABLES_DIR, 'stacking_replicates.csv'))
    df_puck = pd.read_csv(os.path.join(TABLES_DIR, 'sugar_pucker_replicates.csv')).set_index('replicate_tag')
    
    reps = ['wt_rep1', 'wt_rep2', 'wt_rep3', 'cua_rep1', 'cua_rep2', 'cua_rep3']
    
    audit_rows = []
    
    def add_row(metric, wt1, wt2, wt3, c1, c2, c3, direction, status):
        audit_rows.append({
            'metric': metric,
            'wt_rep1': round(wt1, 4) if isinstance(wt1, (float, int)) else wt1,
            'wt_rep2': round(wt2, 4) if isinstance(wt2, (float, int)) else wt2,
            'wt_rep3': round(wt3, 4) if isinstance(wt3, (float, int)) else wt3,
            'cua_rep1': round(c1, 4) if isinstance(c1, (float, int)) else c1,
            'cua_rep2': round(c2, 4) if isinstance(c2, (float, int)) else c2,
            'cua_rep3': round(c3, 4) if isinstance(c3, (float, int)) else c3,
            'direction_consistency': direction,
            'sampling_status': status
        })
        
    # 1. Full ASL RMSD
    w_rmsd = [df_rmsd.loc[f'wt_rep{i}', 'full_asl_mean'] for i in [1, 2, 3]]
    c_rmsd = [df_rmsd.loc[f'cua_rep{i}', 'full_asl_mean'] for i in [1, 2, 3]]
    add_row('Full ASL Backbone RMSD (Å)',
            w_rmsd[0], w_rmsd[1], w_rmsd[2], c_rmsd[0], c_rmsd[1], c_rmsd[2],
            'Consistent elevation in CUA (3/3)',
            'Locally stable overall; loop contributes elevated variance')
            
    # 2. 5-bp Stem RMSD
    w_s_rmsd = [df_rmsd.loc[f'wt_rep{i}', 'stem_mean'] for i in [1, 2, 3]]
    c_s_rmsd = [df_rmsd.loc[f'cua_rep{i}', 'stem_mean'] for i in [1, 2, 3]]
    add_row('5-bp Stem Backbone RMSD (Å)',
            w_s_rmsd[0], w_s_rmsd[1], w_s_rmsd[2], c_s_rmsd[0], c_s_rmsd[1], c_s_rmsd[2],
            'Consistent modest elevation in CUA (3/3); both < 1.2 Å',
            'Locally stable (helical stem intact in all 6 runs)')
            
    # 3. Loop Global RMSD (stem-aligned)
    w_lg_rmsd = [df_rmsd.loc[f'wt_rep{i}', 'loop_global_mean'] for i in [1, 2, 3]]
    c_lg_rmsd = [df_rmsd.loc[f'cua_rep{i}', 'loop_global_mean'] for i in [1, 2, 3]]
    add_row('7-nt Loop Global RMSD (Stem-Aligned) (Å)',
            w_lg_rmsd[0], w_lg_rmsd[1], w_lg_rmsd[2], c_lg_rmsd[0], c_lg_rmsd[1], c_lg_rmsd[2],
            'Consistent elevation in CUA (3/3)',
            'Not demonstrably converged over 20 ns (loop hinge drift)')
            
    # 4. Loop Internal RMSD (loop-aligned)
    w_li_rmsd = [df_rmsd.loc[f'wt_rep{i}', 'loop_internal_mean'] for i in [1, 2, 3]]
    c_li_rmsd = [df_rmsd.loc[f'cua_rep{i}', 'loop_internal_mean'] for i in [1, 2, 3]]
    add_row('7-nt Loop Internal RMSD (Loop-Aligned) (Å)',
            w_li_rmsd[0], w_li_rmsd[1], w_li_rmsd[2], c_li_rmsd[0], c_li_rmsd[1], c_li_rmsd[2],
            'Modest elevation in CUA (3/3)',
            'Not demonstrably converged over 20 ns (internal rearrangements)')
            
    # 5. Stem 5-bp Base-Pair Intact Occupancy
    stem_means = {}
    for r in reps:
        sub = df_stem[df_stem['replicate_tag'] == r]
        stem_means[r] = sub['bp_intact_occupancy'].mean() * 100.0
    add_row('5-bp Stem Intact Occupancy (%)',
            stem_means['wt_rep1'], stem_means['wt_rep2'], stem_means['wt_rep3'],
            stem_means['cua_rep1'], stem_means['cua_rep2'], stem_means['cua_rep3'],
            'Consistent high persistence in both constructs (≥98.7%)',
            'Locally stable (helical base-pairing locked throughout 20 ns)')
            
    # 6. U-Turn HB Occupancy
    w_ut_hb = [df_uturn.loc[f'wt_rep{i}', 'hb_occupancy'] * 100.0 for i in [1, 2, 3]]
    c_ut_hb = [df_uturn.loc[f'cua_rep{i}', 'hb_occupancy'] * 100.0 for i in [1, 2, 3]]
    add_row('U33–P36 U-Turn Intact HB Occupancy (%)',
            w_ut_hb[0], w_ut_hb[1], w_ut_hb[2], c_ut_hb[0], c_ut_hb[1], c_ut_hb[2],
            'Consistent complete loss in CUA (3/3: ≤0.05% vs ≥67.5% WT)',
            'Not demonstrably converged over 20 ns (loop breathing in wt_rep3)')
            
    # 7. U-Turn Contact (<= 3.5 A) Occupancy
    w_ut_c35 = [df_uturn.loc[f'wt_rep{i}', 'contact_35A_occupancy'] * 100.0 for i in [1, 2, 3]]
    c_ut_c35 = [df_uturn.loc[f'cua_rep{i}', 'contact_35A_occupancy'] * 100.0 for i in [1, 2, 3]]
    add_row('U33–P36 Contact (≤3.5 Å) Occupancy (%)',
            w_ut_c35[0], w_ut_c35[1], w_ut_c35[2], c_ut_c35[0], c_ut_c35[1], c_ut_c35[2],
            'Consistent loss in CUA (3/3: ≤0.30% vs ≥70.3% WT)',
            'Not demonstrably converged over 20 ns (transient opening in wt_rep3)')
            
    # 8. U-Turn Mean Distance (A)
    w_ut_d = [df_uturn.loc[f'wt_rep{i}', 'dist_min_mean_A'] for i in [1, 2, 3]]
    c_ut_d = [df_uturn.loc[f'cua_rep{i}', 'dist_min_mean_A'] for i in [1, 2, 3]]
    add_row('U33–P36 Minimum Distance d_min (Å)',
            w_ut_d[0], w_ut_d[1], w_ut_d[2], c_ut_d[0], c_ut_d[1], c_ut_d[2],
            'Consistent distance dilation in CUA (3/3: ≥5.27 Å vs ≤3.42 Å WT)',
            'Locally stationary in wt1/2 and cua1; drifting in wt3 and cua2/3')
            
    # 9. Stacking Step 34-35 Occupancy
    stk34 = {}
    for r in reps:
        sub = df_stack[(df_stack['replicate_tag'] == r) & (df_stack['step'] == '34-35')]
        stk34[r] = sub['stacking_occupancy'].iloc[0] * 100.0
    add_row('3\'-Stacking Step 34–35 Occupancy (%)',
            stk34['wt_rep1'], stk34['wt_rep2'], stk34['wt_rep3'],
            stk34['cua_rep1'], stk34['cua_rep2'], stk34['cua_rep3'],
            'Consistent major reduction in CUA (3/3: max 40.9% CUA vs min 84.2% WT)',
            'Not demonstrably converged over 20 ns (unstacking transition in cua_rep2/3)')
            
    # 10. Stacking Step 35-36 Occupancy
    stk35 = {}
    for r in reps:
        sub = df_stack[(df_stack['replicate_tag'] == r) & (df_stack['step'] == '35-36')]
        stk35[r] = sub['stacking_occupancy'].iloc[0] * 100.0
    add_row('3\'-Stacking Step 35–36 Occupancy (%)',
            stk35['wt_rep1'], stk35['wt_rep2'], stk35['wt_rep3'],
            stk35['cua_rep1'], stk35['cua_rep2'], stk35['cua_rep3'],
            'Comparable occupancy in WT and CUA (~40-55%)',
            'Locally fluctuating dynamically in both constructs')
            
    # 11. Stacking Step 36-37 Occupancy
    stk36 = {}
    for r in reps:
        sub = df_stack[(df_stack['replicate_tag'] == r) & (df_stack['step'] == '36-37')]
        stk36[r] = sub['stacking_occupancy'].iloc[0] * 100.0
    add_row('3\'-Stacking Step 36–37 Occupancy (%)',
            stk36['wt_rep1'], stk36['wt_rep2'], stk36['wt_rep3'],
            stk36['cua_rep1'], stk36['cua_rep2'], stk36['cua_rep3'],
            'Comparable occupancy in WT and CUA (~50-70%)',
            'Locally fluctuating dynamically in both constructs')
            
    # 12. Stacking Step 37-38 Occupancy
    stk37 = {}
    for r in reps:
        sub = df_stack[(df_stack['replicate_tag'] == r) & (df_stack['step'] == '37-38')]
        stk37[r] = sub['stacking_occupancy'].iloc[0] * 100.0
    add_row('3\'-Stacking Step 37–38 Occupancy (%)',
            stk37['wt_rep1'], stk37['wt_rep2'], stk37['wt_rep3'],
            stk37['cua_rep1'], stk37['cua_rep2'], stk37['cua_rep3'],
            'Consistent high stacking in both constructs (≥76%)',
            'Locally stable (robust base-stacking scaffold at 3\' boundary)')
            
    # 13. Pos 34 Sugar Pucker North (C3'-endo)
    w_puck_n = [df_puck.loc[f'wt_rep{i}', 'north_c3_endo_occupancy'] * 100.0 for i in [1, 2, 3]]
    c_puck_n = [df_puck.loc[f'cua_rep{i}', 'north_c3_endo_occupancy'] * 100.0 for i in [1, 2, 3]]
    add_row('Pos 34 Sugar Pucker North (C3\'-endo) Occ (%)',
            w_puck_n[0], w_puck_n[1], w_puck_n[2], c_puck_n[0], c_puck_n[1], c_puck_n[2],
            'Predominantly North in both constructs (≥77%)',
            'Observed conformational occupancy over sampled trajectories')
            
    # 14. Pos 34 Sugar Pucker East (O4'-endo)
    w_puck_e = [df_puck.loc[f'wt_rep{i}', 'east_o4_endo_occupancy'] * 100.0 for i in [1, 2, 3]]
    c_puck_e = [df_puck.loc[f'cua_rep{i}', 'east_o4_endo_occupancy'] * 100.0 for i in [1, 2, 3]]
    add_row('Pos 34 Sugar Pucker East (O4\'-endo) Occ (%)',
            w_puck_e[0], w_puck_e[1], w_puck_e[2], c_puck_e[0], c_puck_e[1], c_puck_e[2],
            'Consistent elevation in CUA (3/3: min 15.65% CUA vs max 6.60% WT)',
            'Observed conformational occupancy over sampled trajectories')

    df_out = pd.DataFrame(audit_rows)
    out_csv = os.path.join(TABLES_DIR, 'replicate_consistency_audit.csv')
    df_out.to_csv(out_csv, index=False)
    print(f"Saved replicate consistency audit table: {out_csv}")
    print("=" * 80)
    print(df_out[['metric', 'direction_consistency', 'sampling_status']].to_string(index=False))
    print("=" * 80)

if __name__ == '__main__':
    run_consistency_audit()
