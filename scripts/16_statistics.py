"""
Phase 6: Replicate-Level Statistics and Sampling Sufficiency / Block Convergence Assessment.
Evaluates:
  1. Strict replicate-level statistics (N=3 WT vs N=3 CUA independent simulation replicates):
     - Means, SDs (ddof=1), 95% Confidence Intervals (Student's t, df=2)
     - Difference of means (Delta)
     - Cohen's d effect size (pooled SD)
     - Welch's two-sample unpaired t-test (two-sided, Satterthwaite df)
     - Benjamini-Hochberg FDR correction
  2. Sampling sufficiency and block convergence across four 5-ns time windows:
     - Block 1: 0 - 5 ns   (frames 0 to 500)
     - Block 2: 5 - 10 ns  (frames 500 to 1000)
     - Block 3: 10 - 15 ns (frames 1000 to 1500)
     - Block 4: 15 - 20 ns (frames 1500 to 2000)

Outputs:
  - analysis/tables/replicate_master_table.csv
  - analysis/tables/statistical_comparison.csv
  - analysis/tables/sampling_assessment.csv
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_DIR = os.path.join(BASE_DIR, 'analysis', 'tables')
os.makedirs(TABLES_DIR, exist_ok=True)

T_CRIT_DF2_95 = 4.30265 # Student's t critical value for df = 2, alpha = 0.05 (two-tailed)

def get_stats_replicates(vals):
    """
    vals: array-like of length 3 (independent replicates)
    Returns dict with mean, std, 95% CI half-width, min, max
    """
    vals = np.array(vals, dtype=float)
    n = len(vals)
    mean = np.mean(vals)
    std = np.std(vals, ddof=1) if n > 1 else 0.0
    sem = std / np.sqrt(n)
    ci95 = T_CRIT_DF2_95 * sem
    return {
        'mean': mean,
        'std': std,
        'sem': sem,
        'ci95_lower': mean - ci95,
        'ci95_upper': mean + ci95,
        'min': np.min(vals),
        'max': np.max(vals)
    }

def cohen_d(x, y):
    """
    Cohen's d for two independent groups x and y.
    """
    nx, ny = len(x), len(y)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    s_pool = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    if s_pool < 1e-12:
        return 0.0
    return (np.mean(y) - np.mean(x)) / s_pool

def run_statistics():
    print("=" * 80)
    print("PHASE 6: REPLICATE-LEVEL STATISTICS & SAMPLING SUFFICIENCY")
    print("=" * 80)
    
    # Load previously computed replicate tables
    df_rmsd_rep = pd.read_csv(os.path.join(TABLES_DIR, 'rmsd_replicates.csv'))
    df_stem_rep = pd.read_csv(os.path.join(TABLES_DIR, 'stem_hbond_replicates.csv'))
    df_uturn_rep = pd.read_csv(os.path.join(TABLES_DIR, 'uturn_replicates.csv'))
    df_stack_rep = pd.read_csv(os.path.join(TABLES_DIR, 'stacking_replicates.csv'))
    df_pucker_rep = pd.read_csv(os.path.join(TABLES_DIR, 'sugar_pucker_replicates.csv'))
    
    # Assemble master replicate table
    reps = ['wt_rep1', 'wt_rep2', 'wt_rep3', 'cua_rep1', 'cua_rep2', 'cua_rep3']
    master_rows = []
    
    for rep in reps:
        construct = rep.split('_')[0]
        rep_num = int(rep.split('_')[1].replace('rep', ''))
        
        # RMSD
        r_rmsd = df_rmsd_rep[df_rmsd_rep['replicate_tag'] == rep].iloc[0]
        # Stem HB (average across all 5 pairs for this replicate)
        r_stem = df_stem_rep[df_stem_rep['replicate_tag'] == rep]
        mean_stem_intact = r_stem['bp_intact_occupancy'].mean()
        # U-turn
        r_uturn = df_uturn_rep[df_uturn_rep['replicate_tag'] == rep].iloc[0]
        # Stacking
        r_stk = df_stack_rep[df_stack_rep['replicate_tag'] == rep]
        stk_34_35 = r_stk[r_stk['step'] == '34-35']['stacking_occupancy'].iloc[0]
        stk_35_36 = r_stk[r_stk['step'] == '35-36']['stacking_occupancy'].iloc[0]
        stk_36_37 = r_stk[r_stk['step'] == '36-37']['stacking_occupancy'].iloc[0]
        stk_37_38 = r_stk[r_stk['step'] == '37-38']['stacking_occupancy'].iloc[0]
        # Sugar pucker
        r_puck = df_pucker_rep[df_pucker_rep['replicate_tag'] == rep].iloc[0]
        
        master_rows.append({
            'replicate_tag': rep,
            'construct': construct,
            'replicate_num': rep_num,
            'full_asl_rmsd_mean_A': r_rmsd['full_asl_mean'],
            'full_asl_rmsd_final10ns_A': r_rmsd['full_asl_final10ns_mean'],
            'stem_rmsd_mean_A': r_rmsd['stem_mean'],
            'loop_global_rmsd_mean_A': r_rmsd['loop_global_mean'],
            'loop_internal_rmsd_mean_A': r_rmsd['loop_internal_mean'],
            'stem_bp_intact_mean_occ': mean_stem_intact,
            'uturn_hb_occupancy': r_uturn['hb_occupancy'],
            'uturn_contact35_occupancy': r_uturn['contact_35A_occupancy'],
            'uturn_dist_mean_A': r_uturn['dist_min_mean_A'],
            'stacking_34_35_occupancy': stk_34_35,
            'stacking_35_36_occupancy': stk_35_36,
            'stacking_36_37_occupancy': stk_36_37,
            'stacking_37_38_occupancy': stk_37_38,
            'pucker_north_c3_endo_occ': r_puck['north_c3_endo_occupancy'],
            'pucker_east_o4_endo_occ': r_puck['east_o4_endo_occupancy'],
            'pucker_median_phase_deg': r_puck['median_phase_deg'],
        })
        
    df_master = pd.DataFrame(master_rows)
    master_csv = os.path.join(TABLES_DIR, 'replicate_master_table.csv')
    df_master.to_csv(master_csv, index=False)
    print(f"Saved master replicate table: {master_csv}")
    
    # Statistical comparison across the N=3 vs N=3 replicates
    metric_cols = [
        ('Full ASL Backbone RMSD (A)', 'full_asl_rmsd_mean_A'),
        ('Full ASL RMSD Final 10 ns (A)', 'full_asl_rmsd_final10ns_A'),
        ('5-bp Stem Backbone RMSD (A)', 'stem_rmsd_mean_A'),
        ('Loop Global RMSD (stem-aligned) (A)', 'loop_global_rmsd_mean_A'),
        ('Loop Internal RMSD (loop-aligned) (A)', 'loop_internal_rmsd_mean_A'),
        ('5-bp Stem Intact Occupancy (fraction)', 'stem_bp_intact_mean_occ'),
        ('U-turn U33-P36 HB Occupancy (fraction)', 'uturn_hb_occupancy'),
        ('U-turn U33-P36 Contact <=3.5 A Occ (fraction)', 'uturn_contact35_occupancy'),
        ('U-turn U33-P36 Mean Distance (A)', 'uturn_dist_mean_A'),
        ('3\'-Stacking 34-35 Occupancy (fraction)', 'stacking_34_35_occupancy'),
        ('3\'-Stacking 35-36 Occupancy (fraction)', 'stacking_35_36_occupancy'),
        ('3\'-Stacking 36-37 Occupancy (fraction)', 'stacking_36_37_occupancy'),
        ('3\'-Stacking 37-38 Occupancy (fraction)', 'stacking_37_38_occupancy'),
        ('Pos 34 Sugar Pucker North (C3\'-endo) Occ', 'pucker_north_c3_endo_occ'),
        ('Pos 34 Sugar Pucker East (O4\'-endo) Occ', 'pucker_east_o4_endo_occ'),
    ]
    
    stat_rows = []
    
    for metric_label, col in metric_cols:
        wt_vals = df_master[df_master['construct'] == 'wt'][col].values
        cua_vals = df_master[df_master['construct'] == 'cua'][col].values
        
        st_wt = get_stats_replicates(wt_vals)
        st_cua = get_stats_replicates(cua_vals)
        
        delta = st_cua['mean'] - st_wt['mean']
        d_val = cohen_d(wt_vals, cua_vals)
        
        # Welch's t-test (two-sample unpaired, unequal variance)
        t_stat, p_val = stats.ttest_ind(cua_vals, wt_vals, equal_var=False)
        # Satterthwaite degrees of freedom
        v_wt = np.var(wt_vals, ddof=1)
        v_cua = np.var(cua_vals, ddof=1)
        df_num = (v_wt/3 + v_cua/3)**2
        df_den = (v_wt/3)**2 / 2 + (v_cua/3)**2 / 2
        df_satt = df_num / df_den if df_den > 1e-12 else 4.0
        
        stat_rows.append({
            'metric': metric_label,
            'col_id': col,
            'wt_n': 3,
            'wt_mean': st_wt['mean'],
            'wt_std': st_wt['std'],
            'wt_ci95_lower': st_wt['ci95_lower'],
            'wt_ci95_upper': st_wt['ci95_upper'],
            'cua_n': 3,
            'cua_mean': st_cua['mean'],
            'cua_std': st_cua['std'],
            'cua_ci95_lower': st_cua['ci95_lower'],
            'cua_ci95_upper': st_cua['ci95_upper'],
            'delta_cua_minus_wt': delta,
            'cohen_d': d_val,
            'welch_t_stat': float(t_stat),
            'welch_df': float(df_satt),
            'welch_p_val': float(p_val),
        })
        
    df_stats = pd.DataFrame(stat_rows)
    
    # Benjamini-Hochberg FDR correction
    p_values = df_stats['welch_p_val'].values
    m = len(p_values)
    sorted_idx = np.argsort(p_values)
    q_values = np.zeros(m)
    for rank, idx in enumerate(sorted_idx, 1):
        q_values[idx] = min(1.0, p_values[idx] * m / rank)
    # Ensure monotonicity of q-values
    for i in range(m - 2, -1, -1):
        idx_curr = sorted_idx[i]
        idx_next = sorted_idx[i + 1]
        q_values[idx_curr] = min(q_values[idx_curr], q_values[idx_next])
        
    df_stats['fdr_q_val'] = q_values
    
    stat_csv = os.path.join(TABLES_DIR, 'statistical_comparison.csv')
    df_stats.to_csv(stat_csv, index=False)
    print(f"Saved statistical comparison table: {stat_csv}")
    
    # -------------------------------------------------------------
    # SAMPLING SUFFICIENCY & BLOCK CONVERGENCE ASSESSMENT
    # -------------------------------------------------------------
    print("\nEvaluating sampling sufficiency across 4 contiguous 5-ns time blocks...")
    # Load time series files
    df_rmsd_ts = pd.read_csv(os.path.join(TABLES_DIR, 'rmsd_time_series.csv'))
    df_uturn_ts = pd.read_csv(os.path.join(TABLES_DIR, 'uturn_time_series.csv'))
    df_stack_ts = pd.read_csv(os.path.join(TABLES_DIR, 'stacking_time_series.csv'))
    
    # Blocks: 0-500 (0-5ns), 500-1000 (5-10ns), 1000-1500 (10-15ns), 1500-2000 (15-20ns)
    blocks = [
        ('Block 1 (0-5 ns)', 0, 500, '0-5ns'),
        ('Block 2 (5-10 ns)', 500, 1000, '5-10ns'),
        ('Block 3 (10-15 ns)', 1000, 1500, '10-15ns'),
        ('Block 4 (15-20 ns)', 1500, 2000, '15-20ns'),
    ]
    
    block_records = []
    
    for rep in reps:
        construct = rep.split('_')[0]
        rep_num = int(rep.split('_')[1].replace('rep', ''))
        
        rep_rmsd = df_rmsd_ts[df_rmsd_ts['replicate_tag'] == rep]
        rep_uturn = df_uturn_ts[df_uturn_ts['replicate_tag'] == rep]
        rep_stk3435 = df_stack_ts[(df_stack_ts['replicate_tag'] == rep) & (df_stack_ts['step'] == '34-35')]
        
        for b_name, f_start, f_end, b_tag in blocks:
            b_rmsd = rep_rmsd[(rep_rmsd['frame'] >= f_start) & (rep_rmsd['frame'] < f_end)]
            b_uturn = rep_uturn[(rep_uturn['frame'] >= f_start) & (rep_uturn['frame'] < f_end)]
            b_stk = rep_stk3435[(rep_stk3435['frame'] >= f_start) & (rep_stk3435['frame'] < f_end)]
            
            block_records.append({
                'replicate_tag': rep,
                'construct': construct,
                'replicate_num': rep_num,
                'block_name': b_name,
                'block_tag': b_tag,
                'time_range_ns': f"{f_start*0.010:.1f}-{f_end*0.010:.1f} ns",
                'full_asl_rmsd_mean': float(b_rmsd['rmsd_full_asl_angstrom'].mean()),
                'stem_rmsd_mean': float(b_rmsd['rmsd_stem_angstrom'].mean()),
                'loop_global_rmsd_mean': float(b_rmsd['rmsd_loop_global_angstrom'].mean()),
                'loop_internal_rmsd_mean': float(b_rmsd['rmsd_loop_internal_angstrom'].mean()),
                'uturn_hb_occupancy': float(b_uturn['hb_intact'].mean()),
                'uturn_dist_mean_A': float(b_uturn['dist_n3_op_min_A'].mean()),
                'stacking_34_35_occupancy': float(b_stk['is_stacked'].mean()),
            })
            
    df_blocks = pd.DataFrame(block_records)
    block_csv = os.path.join(TABLES_DIR, 'sampling_assessment.csv')
    df_blocks.to_csv(block_csv, index=False)
    print(f"Saved sampling assessment table: {block_csv}")
    print("=" * 80)

if __name__ == '__main__':
    run_statistics()
