"""
Phase 6: Publication-Quality Figure Generation.
Generates:
  - figures/figure1_rmsd.png: Full-ASL, stem, and loop RMSD time series and distributions.
  - figures/figure2_rmsf.png: Per-residue RMSF profile across all 17 ASL nucleotides.
  - figures/figure3_stem_hbond.png: Stem base-pair hydrogen-bond occupancy.
  - figures/figure4_uturn_diagnostic.png: U-turn diagnostic (U33 N3 ... P36) time series and distance distribution.
  - figures/figure5_stacking.png: 3'-anticodon base stacking geometry and occupancy.
  - figures/figure6_sugar_pucker.png: Nucleotide 34 ribose sugar pucker pseudorotation and state occupancy.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLES_DIR = os.path.join(BASE_DIR, 'analysis', 'tables')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# Publication styling
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titlesize'] = 13
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.0

# Color palette
WT_COLORS = ['#1f77b4', '#4a90e2', '#7fb3d5'] # Blue shades
CUA_COLORS = ['#d95f02', '#e67e22', '#f39c12'] # Orange/Amber shades
WT_MAIN = '#1f77b4'
CUA_MAIN = '#d95f02'

def plot_figure1_rmsd():
    print("Generating Figure 1: Backbone RMSD...")
    df_ts = pd.read_csv(os.path.join(TABLES_DIR, 'rmsd_time_series.csv'))
    df_rep = pd.read_csv(os.path.join(TABLES_DIR, 'rmsd_replicates.csv'))
    
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    
    rep_styles = {
        'wt_rep1': (WT_COLORS[0], '-'),
        'wt_rep2': (WT_COLORS[1], '--'),
        'wt_rep3': (WT_COLORS[2], '-.'),
        'cua_rep1': (CUA_COLORS[0], '-'),
        'cua_rep2': (CUA_COLORS[1], '--'),
        'cua_rep3': (CUA_COLORS[2], '-.'),
    }
    
    # Downsample time-series slightly for smooth plotting (stride 2)
    df_plot = df_ts.iloc[::2]
    
    # (A) Full ASL RMSD
    ax = axes[0, 0]
    for rep, (col, ls) in rep_styles.items():
        sub = df_plot[df_plot['replicate_tag'] == rep]
        lbl = f"WT rep {rep[-1]}" if "wt" in rep else f"CUA rep {rep[-1]}"
        ax.plot(sub['time_ns'], sub['rmsd_full_asl_angstrom'], label=lbl, color=col, linestyle=ls, alpha=0.85, linewidth=1.2)
    ax.set_title('A. Full 17-nt ASL Backbone RMSD', fontweight='bold', loc='left')
    ax.set_xlabel('Simulation Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_ylim(0, 4.5)
    ax.legend(loc='upper right', ncol=2, frameon=True, framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # (B) Stem RMSD
    ax = axes[0, 1]
    for rep, (col, ls) in rep_styles.items():
        sub = df_plot[df_plot['replicate_tag'] == rep]
        ax.plot(sub['time_ns'], sub['rmsd_stem_angstrom'], color=col, linestyle=ls, alpha=0.85, linewidth=1.2)
    ax.set_title('B. 5-bp Stem Backbone RMSD', fontweight='bold', loc='left')
    ax.set_xlabel('Simulation Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_ylim(0, 2.5)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # (C) Loop Global RMSD (stem-aligned)
    ax = axes[1, 0]
    for rep, (col, ls) in rep_styles.items():
        sub = df_plot[df_plot['replicate_tag'] == rep]
        ax.plot(sub['time_ns'], sub['rmsd_loop_global_angstrom'], color=col, linestyle=ls, alpha=0.85, linewidth=1.2)
    ax.set_title('C. 7-nt Loop Global RMSD (Stem-Aligned)', fontweight='bold', loc='left')
    ax.set_xlabel('Simulation Time (ns)')
    ax.set_ylabel('Displacement (Å)')
    ax.set_ylim(0, 8.5)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # (D) Replicate Summary Comparison
    ax = axes[1, 1]
    metrics = [
        ('Full ASL', 'full_asl_mean'),
        ('Stem', 'stem_mean'),
        ('Loop Global', 'loop_global_mean'),
        ('Loop Internal', 'loop_internal_mean')
    ]
    x_positions = np.arange(len(metrics))
    width = 0.35
    
    wt_means = [df_rep[df_rep['construct'] == 'wt'][m[1]].mean() for m in metrics]
    wt_stds = [df_rep[df_rep['construct'] == 'wt'][m[1]].std(ddof=1) for m in metrics]
    cua_means = [df_rep[df_rep['construct'] == 'cua'][m[1]].mean() for m in metrics]
    cua_stds = [df_rep[df_rep['construct'] == 'cua'][m[1]].std(ddof=1) for m in metrics]
    
    rects1 = ax.bar(x_positions - width/2, wt_means, width, yerr=wt_stds, label='WT (N=3)',
                    color=WT_MAIN, capsize=4, edgecolor='black', alpha=0.85)
    rects2 = ax.bar(x_positions + width/2, cua_means, width, yerr=cua_stds, label='CUA (N=3)',
                    color=CUA_MAIN, capsize=4, edgecolor='black', alpha=0.85)
    
    # Overlay individual replicate dots
    for i, m in enumerate(metrics):
        wt_pts = df_rep[df_rep['construct'] == 'wt'][m[1]].values
        cua_pts = df_rep[df_rep['construct'] == 'cua'][m[1]].values
        ax.scatter([x_positions[i] - width/2]*3, wt_pts, color='navy', zorder=5, s=25, edgecolor='white')
        ax.scatter([x_positions[i] + width/2]*3, cua_pts, color='darkred', zorder=5, s=25, edgecolor='white')
        
    ax.set_title('D. Replicate Mean Backbone RMSD (Mean ± SD, N=3)', fontweight='bold', loc='left')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([m[0] for m in metrics])
    ax.set_ylabel('Backbone RMSD (Å)')
    ax.set_ylim(0, 6.5)
    ax.legend(loc='upper left', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.6, axis='y')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'figure1_rmsd.png')
    plt.savefig(out_path)
    plt.close()
    print(f"Saved: {out_path}")

def plot_figure2_rmsf():
    print("Generating Figure 2: Per-Residue RMSF...")
    df_summary = pd.read_csv(os.path.join(TABLES_DIR, 'rmsf_summary.csv'))
    df_rep = pd.read_csv(os.path.join(TABLES_DIR, 'rmsf_replicates.csv'))
    
    fig, ax = plt.subplots(figsize=(10, 5.5))
    
    wt_sum = df_summary[df_summary['construct'] == 'wt'].sort_values('trna_position')
    cua_sum = df_summary[df_summary['construct'] == 'cua'].sort_values('trna_position')
    
    pos = wt_sum['trna_position'].values
    x = np.arange(len(pos))
    
    # Shaded structural regions
    # Stem: 27-31 (indices 0..4) and 39-43 (indices 12..16)
    # Loop: 32-38 (indices 5..11)
    # Anticodon: 34-36 (indices 7..9)
    ax.axvspan(-0.5, 4.5, color='gray', alpha=0.15, label='5\'-Stem')
    ax.axvspan(4.5, 11.5, color='gold', alpha=0.15, label='7-nt Loop')
    ax.axvspan(6.5, 9.5, color='orange', alpha=0.25, label='Anticodon (34-36)')
    ax.axvspan(11.5, 16.5, color='gray', alpha=0.15, label='3\'-Stem')
    
    # Plot WT mean ± SD
    ax.plot(x, wt_sum['backbone_rmsf_mean'], marker='o', color=WT_MAIN, linewidth=2, label='WT Backbone (Mean ± SD, N=3)')
    ax.fill_between(x, wt_sum['backbone_rmsf_mean'] - wt_sum['backbone_rmsf_std'],
                    wt_sum['backbone_rmsf_mean'] + wt_sum['backbone_rmsf_std'],
                    color=WT_MAIN, alpha=0.25)
    
    # Plot CUA mean ± SD
    ax.plot(x, cua_sum['backbone_rmsf_mean'], marker='s', color=CUA_MAIN, linewidth=2, label='CUA Backbone (Mean ± SD, N=3)')
    ax.fill_between(x, cua_sum['backbone_rmsf_mean'] - cua_sum['backbone_rmsf_std'],
                    cua_sum['backbone_rmsf_mean'] + cua_sum['backbone_rmsf_std'],
                    color=CUA_MAIN, alpha=0.25)
    
    # Individual replicate points
    for rep in [1, 2, 3]:
        sub_wt = df_rep[(df_rep['construct'] == 'wt') & (df_rep['replicate'] == rep)].sort_values('trna_position')
        sub_cua = df_rep[(df_rep['construct'] == 'cua') & (df_rep['replicate'] == rep)].sort_values('trna_position')
        ax.scatter(x, sub_wt['backbone_rmsf_angstrom'], color='navy', alpha=0.4, s=20)
        ax.scatter(x, sub_cua['backbone_rmsf_angstrom'], color='darkred', alpha=0.4, s=20)
        
    labels = [f"{w_res}{p}" if w_res == c_res else f"{w_res}/{c_res}{p}"
              for p, w_res, c_res in zip(pos, wt_sum['resname'], cua_sum['resname'])]
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontweight='medium')
    ax.set_xlabel('Residue & tRNA Canonical Position (5\' -> 3\')')
    ax.set_ylabel('Backbone RMSF (Å)')
    ax.set_title('Per-Residue Backbone Fluctuation (RMSF) Aligned to Stem Reference', fontweight='bold', loc='left')
    ax.set_ylim(0.5, 3.0)
    ax.legend(loc='upper right', frameon=True, framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'figure2_rmsf.png')
    plt.savefig(out_path)
    plt.close()
    print(f"Saved: {out_path}")

def plot_figure3_stem_hbond():
    print("Generating Figure 3: Stem Base-Pair H-Bond Occupancy...")
    df_rep = pd.read_csv(os.path.join(TABLES_DIR, 'stem_hbond_replicates.csv'))
    df_sum = pd.read_csv(os.path.join(TABLES_DIR, 'stem_hbond_summary.csv'))
    
    pairs = ['C27:G43', 'C28:G42', 'C29:G41', 'G30:C40', 'C31:G39']
    x = np.arange(len(pairs))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(9, 5.5))
    
    wt_sum = df_sum[df_sum['construct'] == 'wt'].set_index('base_pair').loc[pairs]
    cua_sum = df_sum[df_sum['construct'] == 'cua'].set_index('base_pair').loc[pairs]
    
    ax.bar(x - width/2, wt_sum['bp_intact_mean']*100, width, yerr=wt_sum['bp_intact_std']*100,
           label='WT (Mean ± SD, N=3)', color=WT_MAIN, capsize=4, edgecolor='black', alpha=0.85)
    ax.bar(x + width/2, cua_sum['bp_intact_mean']*100, width, yerr=cua_sum['bp_intact_std']*100,
           label='CUA (Mean ± SD, N=3)', color=CUA_MAIN, capsize=4, edgecolor='black', alpha=0.85)
    
    # Overlay individual replicates
    for i, pair in enumerate(pairs):
        wt_pts = df_rep[(df_rep['construct'] == 'wt') & (df_rep['base_pair'] == pair)]['bp_intact_occupancy'].values * 100
        cua_pts = df_rep[(df_rep['construct'] == 'cua') & (df_rep['base_pair'] == pair)]['bp_intact_occupancy'].values * 100
        ax.scatter([x[i] - width/2]*3, wt_pts, color='navy', zorder=5, s=25, edgecolor='white')
        ax.scatter([x[i] + width/2]*3, cua_pts, color='darkred', zorder=5, s=25, edgecolor='white')
        
    ax.set_xticks(x)
    ax.set_xticklabels(pairs, fontweight='medium')
    ax.set_xlabel('Stem Base Pair (Watson-Crick G-C)')
    ax.set_ylabel('Base-Pair Intact Occupancy (%)')
    ax.set_title('Stem Base-Pair Stability (>= 2 Watson-Crick Hydrogen Bonds Intact)', fontweight='bold', loc='left')
    ax.set_ylim(90, 102)
    ax.legend(loc='lower left', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.6, axis='y')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'figure3_stem_hbond.png')
    plt.savefig(out_path)
    plt.close()
    print(f"Saved: {out_path}")

def plot_figure4_uturn_diagnostic():
    print("Generating Figure 4: Loop U-Turn Diagnostic...")
    df_ts = pd.read_csv(os.path.join(TABLES_DIR, 'uturn_time_series.csv'))
    df_rep = pd.read_csv(os.path.join(TABLES_DIR, 'uturn_replicates.csv'))
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    
    rep_styles = {
        'wt_rep1': (WT_COLORS[0], '-'),
        'wt_rep2': (WT_COLORS[1], '--'),
        'wt_rep3': (WT_COLORS[2], '-.'),
        'cua_rep1': (CUA_COLORS[0], '-'),
        'cua_rep2': (CUA_COLORS[1], '--'),
        'cua_rep3': (CUA_COLORS[2], '-.'),
    }
    
    # (A) Time series of U33(N3) ... P36(OP) distance
    ax = axes[0]
    df_plot = df_ts.iloc[::2]
    for rep, (col, ls) in rep_styles.items():
        sub = df_plot[df_plot['replicate_tag'] == rep]
        lbl = f"WT rep {rep[-1]}" if "wt" in rep else f"CUA rep {rep[-1]}"
        ax.plot(sub['time_ns'], sub['dist_n3_op_min_A'], label=lbl, color=col, linestyle=ls, alpha=0.85, linewidth=1.1)
    ax.axhline(3.5, color='black', linestyle=':', linewidth=1.2, label='Contact Threshold (3.5 Å)')
    ax.set_title('A. U33(N3) ... P36(OP) Distance Time Series', fontweight='bold', loc='left')
    ax.set_xlabel('Simulation Time (ns)')
    ax.set_ylabel('Distance (Å)')
    ax.set_ylim(2.0, 10.0)
    ax.legend(loc='upper right', ncol=2, fontsize=8, framealpha=0.9)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # (B) Distance Distribution (KDE)
    ax = axes[1]
    wt_dists = df_ts[df_ts['construct'] == 'wt']['dist_n3_op_min_A']
    cua_dists = df_ts[df_ts['construct'] == 'cua']['dist_n3_op_min_A']
    
    sns.kdeplot(wt_dists, ax=ax, color=WT_MAIN, fill=True, alpha=0.35, linewidth=2, label='WT (pooled 60 ns)')
    sns.kdeplot(cua_dists, ax=ax, color=CUA_MAIN, fill=True, alpha=0.35, linewidth=2, label='CUA (pooled 60 ns)')
    ax.axvline(3.5, color='black', linestyle=':', linewidth=1.2, label='Threshold (3.5 Å)')
    ax.set_title('B. Contact Distance Distribution', fontweight='bold', loc='left')
    ax.set_xlabel('Minimum U33(N3) ... P36(OP) Distance (Å)')
    ax.set_ylabel('Probability Density')
    ax.set_xlim(2.0, 10.0)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # (C) HB Occupancy by Replicate
    ax = axes[2]
    x_pos = np.arange(2)
    width = 0.45
    
    wt_occ = df_rep[df_rep['construct'] == 'wt']['hb_occupancy'].values * 100
    cua_occ = df_rep[df_rep['construct'] == 'cua']['hb_occupancy'].values * 100
    
    ax.bar(0, np.mean(wt_occ), width, yerr=np.std(wt_occ, ddof=1), color=WT_MAIN, capsize=5, edgecolor='black', alpha=0.85, label='WT (N=3)')
    ax.bar(1, np.mean(cua_occ), width, yerr=np.std(cua_occ, ddof=1), color=CUA_MAIN, capsize=5, edgecolor='black', alpha=0.85, label='CUA (N=3)')
    
    ax.scatter([0]*3, wt_occ, color='navy', zorder=5, s=35, edgecolor='white')
    ax.scatter([1]*3, cua_occ, color='darkred', zorder=5, s=35, edgecolor='white')
    
    ax.set_title('C. Intact U-Turn HB Occupancy', fontweight='bold', loc='left')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(['WT\n(U33...C36)', 'CUA\n(U33...A36)'], fontweight='medium')
    ax.set_ylabel('Hydrogen-Bond Occupancy (%)')
    ax.set_ylim(-5, 115)
    ax.grid(True, linestyle=':', alpha=0.6, axis='y')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'figure4_uturn_diagnostic.png')
    plt.savefig(out_path)
    plt.close()
    print(f"Saved: {out_path}")

def plot_figure5_stacking():
    print("Generating Figure 5: 3'-Anticodon Stacking...")
    df_rep = pd.read_csv(os.path.join(TABLES_DIR, 'stacking_replicates.csv'))
    df_sum = pd.read_csv(os.path.join(TABLES_DIR, 'stacking_summary.csv'))
    
    steps = ['34-35', '35-36', '36-37', '37-38']
    step_labels = ['Pos 34-35\n(Anticodon 1-2)', 'Pos 35-36\n(Anticodon 2-3)', 'Pos 36-37\n(Anticodon 3 - A37)', 'Pos 37-38\n(A37 - C38)']
    x = np.arange(len(steps))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    
    wt_sum = df_sum[df_sum['construct'] == 'wt'].set_index('step').loc[steps]
    cua_sum = df_sum[df_sum['construct'] == 'cua'].set_index('step').loc[steps]
    
    ax.bar(x - width/2, wt_sum['stacking_occupancy_mean']*100, width, yerr=wt_sum['stacking_occupancy_std']*100,
           label='WT (Mean ± SD, N=3)', color=WT_MAIN, capsize=4, edgecolor='black', alpha=0.85)
    ax.bar(x + width/2, cua_sum['stacking_occupancy_mean']*100, width, yerr=cua_sum['stacking_occupancy_std']*100,
           label='CUA (Mean ± SD, N=3)', color=CUA_MAIN, capsize=4, edgecolor='black', alpha=0.85)
    
    for i, step in enumerate(steps):
        wt_pts = df_rep[(df_rep['construct'] == 'wt') & (df_rep['step'] == step)]['stacking_occupancy'].values * 100
        cua_pts = df_rep[(df_rep['construct'] == 'cua') & (df_rep['step'] == step)]['stacking_occupancy'].values * 100
        ax.scatter([x[i] - width/2]*3, wt_pts, color='navy', zorder=5, s=30, edgecolor='white')
        ax.scatter([x[i] + width/2]*3, cua_pts, color='darkred', zorder=5, s=30, edgecolor='white')
        
    ax.set_xticks(x)
    ax.set_xticklabels(step_labels, fontweight='medium')
    ax.set_xlabel('Loop Dinucleotide Step')
    ax.set_ylabel('Stacking Occupancy (%) [R <= 4.5 Å, θ <= 30°]')
    ax.set_title('3\'-Anticodon Loop Base Stacking Geometry & Occupancy', fontweight='bold', loc='left')
    ax.set_ylim(0, 115)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.6, axis='y')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'figure5_stacking.png')
    plt.savefig(out_path)
    plt.close()
    print(f"Saved: {out_path}")

def plot_figure6_sugar_pucker():
    print("Generating Figure 6: Nucleotide 34 Sugar Pucker...")
    df_ts = pd.read_csv(os.path.join(TABLES_DIR, 'sugar_pucker_time_series.csv'))
    df_rep = pd.read_csv(os.path.join(TABLES_DIR, 'sugar_pucker_replicates.csv'))
    df_sum = pd.read_csv(os.path.join(TABLES_DIR, 'sugar_pucker_summary.csv'))
    
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0))
    
    # (A) Pseudorotation phase angle distribution
    ax = axes[0]
    wt_phases = df_ts[df_ts['construct'] == 'wt']['phase_P_deg']
    cua_phases = df_ts[df_ts['construct'] == 'cua']['phase_P_deg']
    
    bins = np.linspace(0, 360, 73) # 5-degree bins
    ax.hist(wt_phases, bins=bins, color=WT_MAIN, alpha=0.55, density=True, label='WT (G34, pooled 60 ns)', edgecolor='none')
    ax.hist(cua_phases, bins=bins, color=CUA_MAIN, alpha=0.55, density=True, label='CUA (C34, pooled 60 ns)', edgecolor='none')
    
    # Annotate pucker zones
    ax.axvspan(0, 36, color='blue', alpha=0.08)
    ax.axvspan(324, 360, color='blue', alpha=0.08)
    ax.text(18, ax.get_ylim()[1]*0.9 if ax.get_ylim()[1] > 0 else 0.015, 'North\n(C3\'-endo)', ha='center', fontsize=8, color='darkblue')
    ax.axvspan(36, 144, color='green', alpha=0.06)
    ax.text(90, 0.012, 'East\n(O4\'-endo)', ha='center', fontsize=8, color='darkgreen')
    ax.axvspan(144, 216, color='red', alpha=0.06)
    ax.text(180, 0.012, 'South\n(C2\'-endo)', ha='center', fontsize=8, color='darkred')
    
    ax.set_title('A. Pseudorotation Phase Angle P Distribution', fontweight='bold', loc='left')
    ax.set_xlabel('Pseudorotation Phase Angle P (degrees)')
    ax.set_ylabel('Probability Density')
    ax.set_xlim(0, 360)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # (B) State occupancy breakdown
    ax = axes[1]
    states = ['North\n(C3\'-endo)', 'East\n(O4\'-endo)', 'South\n(C2\'-endo)', 'West\n(C1\'-exo)']
    x = np.arange(len(states))
    width = 0.35
    
    wt_occ = [
        df_rep[df_rep['construct'] == 'wt']['north_c3_endo_occupancy'].mean() * 100,
        df_rep[df_rep['construct'] == 'wt']['east_o4_endo_occupancy'].mean() * 100,
        df_rep[df_rep['construct'] == 'wt']['south_c2_endo_occupancy'].mean() * 100,
        df_rep[df_rep['construct'] == 'wt']['west_c1_exo_occupancy'].mean() * 100,
    ]
    wt_std = [
        df_rep[df_rep['construct'] == 'wt']['north_c3_endo_occupancy'].std(ddof=1) * 100,
        df_rep[df_rep['construct'] == 'wt']['east_o4_endo_occupancy'].std(ddof=1) * 100,
        df_rep[df_rep['construct'] == 'wt']['south_c2_endo_occupancy'].std(ddof=1) * 100,
        df_rep[df_rep['construct'] == 'wt']['west_c1_exo_occupancy'].std(ddof=1) * 100,
    ]
    
    cua_occ = [
        df_rep[df_rep['construct'] == 'cua']['north_c3_endo_occupancy'].mean() * 100,
        df_rep[df_rep['construct'] == 'cua']['east_o4_endo_occupancy'].mean() * 100,
        df_rep[df_rep['construct'] == 'cua']['south_c2_endo_occupancy'].mean() * 100,
        df_rep[df_rep['construct'] == 'cua']['west_c1_exo_occupancy'].mean() * 100,
    ]
    cua_std = [
        df_rep[df_rep['construct'] == 'cua']['north_c3_endo_occupancy'].std(ddof=1) * 100,
        df_rep[df_rep['construct'] == 'cua']['east_o4_endo_occupancy'].std(ddof=1) * 100,
        df_rep[df_rep['construct'] == 'cua']['south_c2_endo_occupancy'].std(ddof=1) * 100,
        df_rep[df_rep['construct'] == 'cua']['west_c1_exo_occupancy'].std(ddof=1) * 100,
    ]
    
    ax.bar(x - width/2, wt_occ, width, yerr=wt_std, label='WT (G34)', color=WT_MAIN, capsize=4, edgecolor='black', alpha=0.85)
    ax.bar(x + width/2, cua_occ, width, yerr=cua_std, label='CUA (C34)', color=CUA_MAIN, capsize=4, edgecolor='black', alpha=0.85)
    
    ax.set_xticks(x)
    ax.set_xticklabels(states, fontweight='medium')
    ax.set_ylabel('Conformational State Occupancy (%)')
    ax.set_title('B. Ribose Conformational State Occupancy', fontweight='bold', loc='left')
    ax.set_ylim(0, 105)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(True, linestyle=':', alpha=0.6, axis='y')
    
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'figure6_sugar_pucker.png')
    plt.savefig(out_path)
    plt.close()
    print(f"Saved: {out_path}")

def run_all_figures():
    plot_figure1_rmsd()
    plot_figure2_rmsf()
    plot_figure3_stem_hbond()
    plot_figure4_uturn_diagnostic()
    plot_figure5_stacking()
    plot_figure6_sugar_pucker()
    print("=" * 80)
    print("ALL 6 PUBLICATION FIGURES SUCCESSFULLY GENERATED IN figures/")
    print("=" * 80)

if __name__ == '__main__':
    run_all_figures()
