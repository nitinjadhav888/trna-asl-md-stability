#!/usr/bin/env python3
"""
scripts/02_build_asl_structures.py - Deterministic 3D ASL Structural Builder

Phase 3 Objective:
Construct stereochemically valid, clash-free starting 3D structures for:
1. Wild-Type 17-nt ASL:        CCCGCCUGUCACGCGGG (Sprinzl 27-43)
2. CUA Suppressor 17-nt ASL:   CCCGCCUCUAACGCGGG (Sprinzl 27-43)

Methodology:
- Experimental crystallographic template: PDB 6UGG (Chain A, residues 27-43, 1.95 A).
- Ground truth loop identity: 6UGG loop (residues 32-38: CUGUCAC) is 100% sequence-identical
  to human cytoplasmic tRNA-Asp-GTC-2-1.
- Deterministic glycosidic base grafting:
  * WT construct: stem compensatory pair mutation U29C and A41G.
  * Suppressor construct: WT backbone + anticodon mutations G34C and C36A.
- Terminal chemistry: Neutral 5'-OH (no P/OP1/OP2 on residue 27) and neutral 3'-OH (SDR-09).
- Output: structures/wt/initial.pdb and structures/cua/initial.pdb.
"""

import os
import sys
import numpy as np
from Bio.PDB import PDBParser

def extract_base_template(res):
    """Extracts local orthonormal coordinates of a base ring from an experimental residue."""
    resname = res.resname
    is_purine = resname in ['A', 'G']
    n_anchor = 'N9' if is_purine else 'N1'
    c_dir = 'C4' if is_purine else 'C2'
    
    if is_purine:
        base_atoms = ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4']
        if resname == 'A':
            base_atoms.append('N6')
        elif resname == 'G':
            base_atoms.extend(['O6', 'N2'])
    else:
        base_atoms = ['N1', 'C2', 'O2', 'N3', 'C4', 'C5', 'C6']
        if resname == 'C':
            base_atoms.append('N4')
        elif resname == 'U':
            base_atoms.append('O4')
            
    coords = {a: res[a].get_coord().copy() for a in base_atoms}
    r_n = coords[n_anchor]
    r_dir = coords[c_dir]
    
    vx = r_dir - r_n
    ex = vx / np.linalg.norm(vx)
    
    if is_purine:
        v_other = coords['C8'] - r_n
    else:
        v_other = coords['C6'] - r_n
    ez = np.cross(ex, v_other)
    ez /= np.linalg.norm(ez)
    ey = np.cross(ez, ex)
    
    local_coords = {}
    for a, pos in coords.items():
        rel = pos - r_n
        local_coords[a] = np.array([np.dot(rel, ex), np.dot(rel, ey), np.dot(rel, ez)])
        
    return resname, n_anchor, c_dir, local_coords

def graft_base(res, target_resname, templates):
    """Rigid-body grafts target_resname onto residue res preserving glycosidic bond and chi orientation."""
    orig_name = res.resname
    orig_is_purine = orig_name in ['A', 'G']
    orig_n_anchor = 'N9' if orig_is_purine else 'N1'
    orig_c_dir = 'C4' if orig_is_purine else 'C2'
    
    c1p = res["C1'"].get_coord()
    r_n_orig = res[orig_n_anchor].get_coord()
    r_dir_orig = res[orig_c_dir].get_coord()
    
    # Standard RNA C1'-N bond length (1.470 A)
    u_glyc = (r_n_orig - c1p) / np.linalg.norm(r_n_orig - c1p)
    d_glyc = 1.470
    r_n_target = c1p + d_glyc * u_glyc
    
    # Plane normal
    if orig_is_purine:
        v_other = res['C8'].get_coord() - r_n_orig
    else:
        v_other = res['C6'].get_coord() - r_n_orig
        
    vx_orig = r_dir_orig - r_n_orig
    ex = vx_orig / np.linalg.norm(vx_orig)
    
    ez = np.cross(ex, v_other)
    ez /= np.linalg.norm(ez)
    ey = np.cross(ez, ex)
    
    t_name, t_n_anchor, t_c_dir, t_local = templates[target_resname]
    new_coords = {}
    for a, lcoord in t_local.items():
        pos = r_n_target + lcoord[0] * ex + lcoord[1] * ey + lcoord[2] * ez
        new_coords[a] = pos
        
    return new_coords

def build_asl_coordinates(chain, mutations, templates):
    """
    Constructs atom coordinate dictionary for residues 27 to 43.
    mutations: dict of {res_id: target_resname}
    """
    asl = {}
    for r in range(27, 44):
        res = chain[r]
        target_name = mutations.get(r, res.resname)
        
        # Sugar-phosphate backbone atoms
        # Residue 27: Neutral 5'-OH terminus (omit P, OP1, OP2)
        # Residues 28-43: Internal phosphodiester linkages (include P, OP1, OP2)
        bb_names = ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
        if r > 27:
            bb_names = ['P', 'OP1', 'OP2'] + bb_names
            
        atom_dict = {}
        for bbn in bb_names:
            atom_dict[bbn] = res[bbn].get_coord().copy()
            
        # Base atoms
        if r in mutations:
            base_coords = graft_base(res, target_name, templates)
            atom_dict.update(base_coords)
        else:
            is_purine = res.resname in ['A', 'G']
            if is_purine:
                b_names = ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4']
                if res.resname == 'A':
                    b_names.append('N6')
                elif res.resname == 'G':
                    b_names.extend(['O6', 'N2'])
            else:
                b_names = ['N1', 'C2', 'O2', 'N3', 'C4', 'C5', 'C6']
                if res.resname == 'C':
                    b_names.append('N4')
                elif res.resname == 'U':
                    b_names.append('O4')
            for bn in b_names:
                atom_dict[bn] = res[bn].get_coord().copy()
                
        asl[r] = (target_name, atom_dict)
    return asl

def write_pdb(asl_dict, out_path, chain_id='A'):
    """Writes standardized PDB file with correct nucleic acid atom order."""
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("REMARK   4 tRNA ASL Structural Stability Study\n")
        f.write(f"REMARK   4 Output File: {os.path.basename(out_path)}\n")
        f.write("REMARK   4 Method: Template-based homology modeling from PDB 6UGG (1.95 A)\n")
        f.write("REMARK   4 Termini: Neutral 5'-OH (res 27) and neutral 3'-OH (res 43)\n")
        atom_serial = 1
        for r_num in sorted(asl_dict.keys()):
            r_name, atoms = asl_dict[r_num]
            p_order = ['P', 'OP1', 'OP2']
            sugar_order = ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"]
            if r_name in ['A', 'G']:
                base_order = ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4']
                if r_name == 'A':
                    base_order.append('N6')
                elif r_name == 'G':
                    base_order.extend(['O6', 'N2'])
            else:
                base_order = ['N1', 'C2', 'O2', 'N3', 'C4', 'C5', 'C6']
                if r_name == 'C':
                    base_order.append('N4')
                elif r_name == 'U':
                    base_order.append('O4')
                    
            atom_list = [a for a in p_order if a in atoms] + [a for a in sugar_order if a in atoms] + [a for a in base_order if a in atoms]
            for aname in atom_list:
                pos = atoms[aname]
                element = aname[0]
                line = f"ATOM  {atom_serial:5d} {aname:<4s} {r_name:>3s} {chain_id}{r_num:4d}    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00 50.00          {element:>2s}\n"
                f.write(line)
                atom_serial += 1
        f.write("TER\nEND\n")
    print(f"Generated PDB: {out_path} ({atom_serial - 1} heavy atoms)")

def main():
    print("=" * 80)
    print("Phase 3: 3D ASL Structural Construction Pipeline")
    print("=" * 80)
    
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    template_pdb = os.path.join(base_dir, 'data', 'raw', '6UGG.pdb')
    if not os.path.exists(template_pdb):
        raise FileNotFoundError(f"Template PDB file not found: {template_pdb}")
        
    print(f"Loading crystallographic template: {template_pdb}")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('6UGG', template_pdb)
    chain = structure[0]['A']
    
    # 1. Extract canonical base templates from 6UGG
    templates = {
        'C': extract_base_template(chain[31]),
        'U': extract_base_template(chain[33]),
        'G': extract_base_template(chain[30]),
        'A': extract_base_template(chain[37])
    }
    print("Extracted canonical base templates (A, C, G, U) from 6UGG.")
    
    # 2. Build Wild-Type ASL
    # WT sequence: CCCGCCUGUCACGCGGG (Sprinzl 27-43)
    # 6UGG has U29 and A41 in stem; mutate to C29 and G41
    wt_mutations = {29: 'C', 41: 'G'}
    wt_asl = build_asl_coordinates(chain, wt_mutations, templates)
    wt_seq = ''.join(wt_asl[r][0] for r in sorted(wt_asl.keys()))
    print(f"\nConstructing WT ASL (17 nt):")
    print(f"  Sequence ($5' \\rightarrow 3'$): {wt_seq}")
    print(f"  Stem 5': {wt_seq[:5]} | Loop: {wt_seq[5:12]} | Stem 3': {wt_seq[12:]}")
    wt_pdb_path = os.path.join(base_dir, 'structures', 'wt', 'initial.pdb')
    write_pdb(wt_asl, wt_pdb_path)
    
    # 3. Build CUA Suppressor ASL
    # Suppressor sequence: CCCGCCUCUAACGCGGG (Sprinzl 27-43)
    # Differs from WT ONLY at G34C and C36A
    cua_mutations = {29: 'C', 41: 'G', 34: 'C', 36: 'A'}
    cua_asl = build_asl_coordinates(chain, cua_mutations, templates)
    cua_seq = ''.join(cua_asl[r][0] for r in sorted(cua_asl.keys()))
    print(f"\nConstructing CUA Suppressor ASL (17 nt):")
    print(f"  Sequence ($5' \\rightarrow 3'$): {cua_seq}")
    print(f"  Stem 5': {cua_seq[:5]} | Loop: {cua_seq[5:12]} | Stem 3': {cua_seq[12:]}")
    cua_pdb_path = os.path.join(base_dir, 'structures', 'cua', 'initial.pdb')
    write_pdb(cua_asl, cua_pdb_path)
    
    # 4. Verify Controlled Comparative Design
    non_mutated = [r for r in range(27, 44) if r not in [34, 36]]
    max_dev = 0.0
    for r in non_mutated:
        atoms_wt = wt_asl[r][1]
        atoms_cua = cua_asl[r][1]
        for a in atoms_wt:
            dev = np.linalg.norm(atoms_wt[a] - atoms_cua[a])
            if dev > max_dev:
                max_dev = dev
    print(f"\nControlled Design Verification:")
    print(f"  Non-mutated residues (15/17) max coordinate divergence: {max_dev:.6f} A")
    if max_dev > 1e-5:
        raise ValueError(f"Divergence detected in non-mutated residues: {max_dev} A!")
    print("  SUCCESS: WT and CUA share strictly identical non-mutated coordinates.")
    print("=" * 80)

if __name__ == '__main__':
    main()

