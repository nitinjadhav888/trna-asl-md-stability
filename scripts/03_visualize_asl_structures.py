#!/usr/bin/env python3
"""
scripts/03_visualize_asl_structures.py - 3D Structural Visualization of Initial ASL Models

Generates:
- figures/wt_asl_initial.png
- figures/cua_asl_initial.png
- figures/asl_structural_comparison.png

Visualizes:
- Stem architecture (5-bp A-form helix, Watson-Crick C:G pairing)
- Anticodon loop architecture (7-nt loop, invariant U33 U-turn, purine A37)
- Sequence orientation (5' -> 3')
- Controlled mutation comparison (G34C, C36A)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser

def parse_pdb_coords(pdb_path):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('ASL', pdb_path)
    chain = struct[0]['A']
    
    residues = {}
    for res in chain:
        rid = res.id[1]
        rname = res.resname.strip()
        atoms = {a.name: a.get_coord() for a in res}
        residues[rid] = (rname, atoms)
    return residues

def render_asl_figure(residues, title, out_path, is_suppressor=False):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    
    fig = plt.figure(figsize=(10, 10), dpi=300)
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')
    
    # 1. Extract C1' and P coordinates for backbone trace
    bb_coords = []
    r_ids = sorted(residues.keys())
    for rid in r_ids:
        rname, atoms = residues[rid]
        # Use C1' as representative sugar anchor
        bb_coords.append(atoms["C1'"])
    bb_coords = np.array(bb_coords)
    
    # Plot smooth backbone spline
    ax.plot(bb_coords[:, 0], bb_coords[:, 1], bb_coords[:, 2], 
            color='#475569', linewidth=2.5, linestyle='-', alpha=0.8, zorder=2)
    
    # Color palette
    # Stem 5': Royal Blue (#2563eb)
    # Stem 3': Cyan (#0284c7)
    # Loop Flanking: Slate (#64748b)
    # U33 Invariant: Emerald Green (#059669)
    # Anticodon: Crimson (#dc2626) if WT, Violet (#7c3aed) if Suppressor
    # A37 Purine: Amber (#d97706)
    
    ac_color = '#7c3aed' if is_suppressor else '#dc2626'
    
    for rid in r_ids:
        rname, atoms = residues[rid]
        c1p = atoms["C1'"]
        
        # Color categorization
        if rid in [27, 28, 29, 30, 31]:
            col = '#2563eb' # Stem 5'
            reg = 'Stem 5\''
        elif rid in [39, 40, 41, 42, 43]:
            col = '#0284c7' # Stem 3'
            reg = 'Stem 3\''
        elif rid == 33:
            col = '#059669' # U-turn
            reg = 'U-turn'
        elif rid in [34, 35, 36]:
            col = ac_color   # Anticodon
            reg = 'Anticodon'
        elif rid == 37:
            col = '#d97706' # Purine 37
            reg = 'Purine 37'
        else:
            col = '#64748b' # Flanking loop
            reg = 'Loop'
            
        # Draw sugar-base connection and base atoms
        is_purine = rname in ['A', 'G']
        n_anc = 'N9' if is_purine else 'N1'
        r_n = atoms[n_anc]
        
        # Line from C1' to base
        ax.plot([c1p[0], r_n[0]], [c1p[1], r_n[1]], [c1p[2], r_n[2]], 
                color=col, linewidth=2.0, alpha=0.9)
        
        # Base ring centroid
        ring_atoms = ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'] if is_purine else ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']
        ring_coords = np.array([atoms[a] for a in ring_atoms if a in atoms])
        centroid = np.mean(ring_coords, axis=0)
        
        # Scatter for C1' and base centroid
        ax.scatter([c1p[0]], [c1p[1]], [c1p[2]], color='#1e293b', s=45, zorder=5)
        ax.scatter([centroid[0]], [centroid[1]], [centroid[2]], color=col, s=80, edgecolors='#0f172a', linewidth=0.8, zorder=6)
        
        # Labeling key positions
        label = f"{rname}{rid}"
        if rid in [27, 31, 33, 34, 35, 36, 37, 39, 43]:
            offset = (centroid - c1p) * 0.4
            pos_txt = centroid + offset
            ax.text(pos_txt[0], pos_txt[1], pos_txt[2], label, 
                    fontsize=8.5, fontweight='bold', color='#0f172a')

    # Draw stem base pairs (Watson-Crick H-bonds)
    stem_pairs = [(27, 43), (28, 42), (29, 41), (30, 40), (31, 39)]
    for r5, r3 in stem_pairs:
        c1p_5 = residues[r5][1]["C1'"]
        c1p_3 = residues[r3][1]["C1'"]
        
        # Draw dashed line between paired C1'
        ax.plot([c1p_5[0], c1p_3[0]], [c1p_5[1], c1p_3[1]], [c1p_5[2], c1p_3[2]], 
                color='#0284c7', linestyle='--', linewidth=1.5, alpha=0.6)
        
    # Draw U-turn diagnostic interaction (U33 N3 to 36 OP2)
    u33_n3 = residues[33][1]['N3']
    p36_op2 = residues[36][1]['OP2']
    d_uturn = np.linalg.norm(u33_n3 - p36_op2)
    ax.plot([u33_n3[0], p36_op2[0]], [u33_n3[1], p36_op2[1]], [u33_n3[2], p36_op2[2]], 
            color='#059669', linestyle=':', linewidth=2.0, alpha=0.8)
    
    # Configure 3D view
    ax.view_init(elev=20, azim=45)
    ax.set_axis_off()
    
    # Header & Legend
    seq_str = ''.join(residues[r][0] for r in r_ids)
    plt.title(f"{title}\nSequence: 5'-{seq_str}-3'\nU-Turn Diagnostic: U33(N3)...36(OP2) = {d_uturn:.2f} Å", 
              fontsize=12, fontweight='bold', pad=20, color='#0f172a')
    
    # Annotations box
    legend_text = (
        "Structural Highlights:\n"
        "• Blue/Cyan: 5-bp Watson-Crick Stem (all C:G)\n"
        "• Green: Invariant U33 (Canonical U-turn)\n"
        f"• {'Violet: CUA Suppressor' if is_suppressor else 'Red: Wild-Type GUC'} Anticodon\n"
        "• Amber: Invariant Purine A37 (3'-Stacking)\n"
        "• Grey: Flanking Pyrimidines (C32, C38)\n"
        "• Termini: Neutral 5'-OH (C27) / 3'-OH (G43)"
    )
    fig.text(0.12, 0.08, legend_text, fontsize=9, 
             bbox=dict(boxstyle='round,pad=0.6', facecolor='#f8fafc', edgecolor='#cbd5e1', alpha=0.95))
    
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Generated visualization: {out_path}")

def render_comparison_figure(wt_res, cua_res, out_path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), subplot_kw={'projection': '3d'}, dpi=300)
    fig.patch.set_facecolor('#ffffff')
    
    for ax, res_dict, is_sup, subtitle, ac_col in [
        (ax1, wt_res, False, "Wild-Type ASL (5'-CCCGCCUGUCACGCGGG-3')", '#dc2626'),
        (ax2, cua_res, True, "CUA Suppressor ASL (5'-CCCGCCUCUAACGCGGG-3')", '#7c3aed')
    ]:
        ax.set_facecolor('#ffffff')
        r_ids = sorted(res_dict.keys())
        bb_coords = np.array([res_dict[r][1]["C1'"] for r in r_ids])
        ax.plot(bb_coords[:, 0], bb_coords[:, 1], bb_coords[:, 2], 
                color='#475569', linewidth=2.5, linestyle='-', alpha=0.8)
        
        for rid in r_ids:
            rname, atoms = res_dict[rid]
            c1p = atoms["C1'"]
            if rid in [27, 28, 29, 30, 31]:
                col = '#2563eb'
            elif rid in [39, 40, 41, 42, 43]:
                col = '#0284c7'
            elif rid == 33:
                col = '#059669'
            elif rid in [34, 35, 36]:
                col = ac_col
            elif rid == 37:
                col = '#d97706'
            else:
                col = '#64748b'
                
            is_purine = rname in ['A', 'G']
            n_anc = 'N9' if is_purine else 'N1'
            r_n = atoms[n_anc]
            ax.plot([c1p[0], r_n[0]], [c1p[1], r_n[1]], [c1p[2], r_n[2]], color=col, linewidth=2.0)
            
            ring_atoms = ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'] if is_purine else ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']
            ring_coords = np.array([atoms[a] for a in ring_atoms if a in atoms])
            centroid = np.mean(ring_coords, axis=0)
            
            ax.scatter([c1p[0]], [c1p[1]], [c1p[2]], color='#1e293b', s=35)
            ax.scatter([centroid[0]], [centroid[1]], [centroid[2]], color=col, s=70, edgecolors='#0f172a', linewidth=0.6)
            
            if rid in [27, 33, 34, 35, 36, 37, 43]:
                offset = (centroid - c1p) * 0.45
                pos_txt = centroid + offset
                ax.text(pos_txt[0], pos_txt[1], pos_txt[2], f"{rname}{rid}", fontsize=8, fontweight='bold', color='#0f172a')
                
        # Stem base pairs
        for r5, r3 in [(27, 43), (28, 42), (29, 41), (30, 40), (31, 39)]:
            c1p_5 = res_dict[r5][1]["C1'"]
            c1p_3 = res_dict[r3][1]["C1'"]
            ax.plot([c1p_5[0], c1p_3[0]], [c1p_5[1], c1p_3[1]], [c1p_5[2], c1p_3[2]], color='#0284c7', linestyle='--', linewidth=1.2, alpha=0.5)
            
        ax.view_init(elev=20, azim=45)
        ax.set_axis_off()
        ax.set_title(subtitle, fontsize=11, fontweight='bold', pad=10)

    fig.suptitle("Initial 3D Conformations: Wild-Type vs. CUA Suppressor 17-nt ASL\n(Template: PDB 6UGG, 1.95 Å | Non-Mutated Backbone RMSD = 0.000 Å)", fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Generated comparison visualization: {out_path}")

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    wt_pdb = os.path.join(base_dir, 'structures', 'wt', 'initial.pdb')
    cua_pdb = os.path.join(base_dir, 'structures', 'cua', 'initial.pdb')
    
    wt_res = parse_pdb_coords(wt_pdb)
    cua_res = parse_pdb_coords(cua_pdb)
    
    render_asl_figure(wt_res, "Wild-Type Human Cytoplasmic tRNA-Asp-GTC-2-1 ASL", 
                      os.path.join(base_dir, 'figures', 'wt_asl_initial.png'), is_suppressor=False)
    render_asl_figure(cua_res, "CUA Amber Suppressor Human tRNA-Asp ASL", 
                      os.path.join(base_dir, 'figures', 'cua_asl_initial.png'), is_suppressor=True)
    render_comparison_figure(wt_res, cua_res, 
                             os.path.join(base_dir, 'figures', 'asl_structural_comparison.png'))
    print("All Phase 3 structural figures successfully generated!")

if __name__ == '__main__':
    main()

