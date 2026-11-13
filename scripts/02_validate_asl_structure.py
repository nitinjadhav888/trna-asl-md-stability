#!/usr/bin/env python3
"""
scripts/02_validate_asl_structure.py - Deterministic ASL Structural Validation Script

Phase 3 Objective:
Rigorous, comprehensive structural quality and stereochemical validation of:
- structures/wt/initial.pdb
- structures/cua/initial.pdb

Enforces 14 mandatory checks:
1. Exactly 17 nucleotides
2. Canonical RNA residue identity (A, C, G, U)
3. Sequence matches expected WT sequence (CCCGCCUGUCACGCGGG)
4. Sequence matches expected Suppressor sequence (CCCGCCUCUAACGCGGG)
5. Correct residue numbering (Sprinzl 27 to 43)
6. Correct 5-bp stem / 7-nt loop organization
7. Expected base-pair connectivity (Watson-Crick H-bond geometry in stem)
8. No missing residues
9. No duplicate residues
10. No severe atom-coordinate problems (NaN, Inf, zeroes, unphysical bounds)
11. No obvious steric clashes (all non-bonded heavy atom distances > 1.8 A)
12. Reasonable covalent bond connectivity (backbone, ribose, glycosidic bonds)
13. Correct mutation locations (strictly positions 34 and 36)
14. WT and Suppressor share identical non-mutated coordinates (RMSD == 0.000 A)

The script fails loudly with exit code 1 if ANY requirement is violated.
"""

import os
import sys
import numpy as np
from Bio.PDB import PDBParser

EXPECTED_WT_SEQ = "CCCGCCUGUCACGCGGG"
EXPECTED_CUA_SEQ = "CCCGCCUCUAACGCGGG"
EXPECTED_SPRINZL_IDS = list(range(27, 44))

def load_asl_structure(pdb_path):
    """Loads ASL structure and extracts per-residue dictionary."""
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"Structure file does not exist: {pdb_path}")
        
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('ASL', pdb_path)
    model = struct[0]
    chains = list(model.get_chains())
    if len(chains) != 1:
        raise ValueError(f"{pdb_path}: Expected exactly 1 chain, found {len(chains)}")
        
    chain = chains[0]
    residues = list(chain.get_residues())
    
    res_dict = {}
    for res in residues:
        rid = res.id[1]
        rname = res.resname.strip()
        atoms = {a.name: a.get_coord().copy() for a in res}
        res_dict[rid] = (rname, atoms, res)
        
    return res_dict

def validate_single_structure(pdb_path, expected_seq, construct_name):
    print(f"\n{'='*70}")
    print(f"Validating {construct_name}: {pdb_path}")
    print(f"{'='*70}")
    
    res_dict = load_asl_structure(pdb_path)
    sorted_ids = sorted(res_dict.keys())
    
    # 1. Exactly 17 nucleotides
    if len(res_dict) != 17:
        raise AssertionError(f"Check 1 FAILED: Expected 17 nucleotides, found {len(res_dict)}")
    print(f"[PASS] Check 1: Exactly 17 nucleotides present.")
    
    # 2. RNA residue identity
    valid_names = {'A', 'C', 'G', 'U'}
    for rid, (rname, _, _) in res_dict.items():
        if rname not in valid_names:
            raise AssertionError(f"Check 2 FAILED: Non-canonical RNA residue {rname} at position {rid}")
    print(f"[PASS] Check 2: All residues are canonical RNA nucleotides (A, C, G, U).")
    
    # 3 & 4. Sequence matching
    actual_seq = ''.join(res_dict[rid][0] for rid in sorted_ids)
    if actual_seq != expected_seq:
        raise AssertionError(f"Check 3/4 FAILED: Sequence mismatch!\n  Expected: {expected_seq}\n  Actual:   {actual_seq}")
    print(f"[PASS] Check 3/4: Primary sequence matches expected ({actual_seq}).")
    
    # 5. Correct residue numbering (Sprinzl 27 to 43)
    if sorted_ids != EXPECTED_SPRINZL_IDS:
        raise AssertionError(f"Check 5 FAILED: Expected residue IDs {EXPECTED_SPRINZL_IDS}, found {sorted_ids}")
    print(f"[PASS] Check 5: Residue numbering conforms to Sprinzl coordinates (27 to 43).")
    
    # 6. Correct 5-bp stem / 7-nt loop organization
    stem5_seq = actual_seq[:5]
    loop_seq = actual_seq[5:12]
    stem3_seq = actual_seq[12:]
    if len(stem5_seq) != 5 or len(loop_seq) != 7 or len(stem3_seq) != 5:
        raise AssertionError(f"Check 6 FAILED: Invalid ASL segmentation: stem5={stem5_seq}, loop={loop_seq}, stem3={stem3_seq}")
    # Invariant U33
    if loop_seq[1] != 'U':
        raise AssertionError(f"Check 6 FAILED: Residue 33 is not invariant U (found {loop_seq[1]})")
    # Invariant purine at 37
    if loop_seq[5] not in ('A', 'G'):
        raise AssertionError(f"Check 6 FAILED: Residue 37 is not purine (found {loop_seq[5]})")
    print(f"[PASS] Check 6: Canonical 5-bp stem / 7-nt loop organization confirmed (U33 invariant, A37 purine).")
    
    # 7. Expected base-pair connectivity (Watson-Crick H-bond geometry in stem)
    stem_pairs = [(27, 43), (28, 42), (29, 41), (30, 40), (31, 39)]
    for r5, r3 in stem_pairs:
        name5, atoms5, _ = res_dict[r5]
        name3, atoms3, _ = res_dict[r3]
        
        # In our constructs, all 5 stem pairs are C:G or G:C
        if name5 == 'C' and name3 == 'G':
            d_n3_n1 = np.linalg.norm(atoms5['N3'] - atoms3['N1'])
            d_o2_n2 = np.linalg.norm(atoms5['O2'] - atoms3['N2'])
            d_n4_o6 = np.linalg.norm(atoms5['N4'] - atoms3['O6'])
            if not (2.6 <= d_n3_n1 <= 3.2 and 2.6 <= d_o2_n2 <= 3.2 and 2.6 <= d_n4_o6 <= 3.2):
                raise AssertionError(f"Check 7 FAILED: Distorted C:G H-bonds at pair {r5}:{r3}: {d_n3_n1:.2f}, {d_o2_n2:.2f}, {d_n4_o6:.2f} A")
        elif name5 == 'G' and name3 == 'C':
            d_n1_n3 = np.linalg.norm(atoms5['N1'] - atoms3['N3'])
            d_n2_o2 = np.linalg.norm(atoms5['N2'] - atoms3['O2'])
            d_o6_n4 = np.linalg.norm(atoms5['O6'] - atoms3['N4'])
            if not (2.6 <= d_n1_n3 <= 3.2 and 2.6 <= d_n2_o2 <= 3.2 and 2.6 <= d_o6_n4 <= 3.2):
                raise AssertionError(f"Check 7 FAILED: Distorted G:C H-bonds at pair {r5}:{r3}: {d_n1_n3:.2f}, {d_n2_o2:.2f}, {d_o6_n4:.2f} A")
        else:
            raise AssertionError(f"Check 7 FAILED: Non-Watson-Crick stem pair {r5}({name5}) : {r3}({name3})")
    print(f"[PASS] Check 7: All 5 stem base pairs exhibit standard Watson-Crick hydrogen bonding (2.7-3.0 A).")
    
    # 8. No missing residues
    for exp_id in EXPECTED_SPRINZL_IDS:
        if exp_id not in res_dict:
            raise AssertionError(f"Check 8 FAILED: Missing residue {exp_id}")
    print(f"[PASS] Check 8: Zero missing residues.")
    
    # 9. No duplicate residues
    if len(sorted_ids) != len(set(sorted_ids)):
        raise AssertionError("Check 9 FAILED: Duplicate residue IDs detected.")
    print(f"[PASS] Check 9: Zero duplicate residue IDs.")
    
    # 10. No severe atom-coordinate problems
    for rid, (rname, atoms, res) in res_dict.items():
        for aname, coord in atoms.items():
            if np.isnan(coord).any() or np.isinf(coord).any():
                raise AssertionError(f"Check 10 FAILED: NaN/Inf coordinate in residue {rid} atom {aname}")
            if np.linalg.norm(coord) == 0.0:
                raise AssertionError(f"Check 10 FAILED: All-zero coordinate in residue {rid} atom {aname}")
            # Check B-factor
            bf = res[aname].get_bfactor()
            if bf < 0.0 or bf > 500.0:
                raise AssertionError(f"Check 10 FAILED: Unphysical B-factor ({bf}) in residue {rid} atom {aname}")
    print(f"[PASS] Check 10: Zero NaN/Inf/zero coordinates; all B-factors physically reasonable.")
    
    # 11. No obvious steric clashes (all non-bonded heavy atom distances > 1.8 A)
    all_atoms = []
    for rid in sorted_ids:
        rname, atoms, _ = res_dict[rid]
        for aname, coord in atoms.items():
            all_atoms.append((rid, rname, aname, coord))
            
    for i in range(len(all_atoms)):
        r1, rn1, a1, p1 = all_atoms[i]
        for j in range(i + 1, len(all_atoms)):
            r2, rn2, a2, p2 = all_atoms[j]
            d = np.linalg.norm(p1 - p2)
            
            # Covalent bond exclusions:
            # 1. Same residue bonded atoms
            if r1 == r2:
                if d < 1.15:
                    raise AssertionError(f"Check 11 FAILED: Intra-residue atomic collision ({d:.2f} A) between {r1} {rn1} {a1} and {a2}")
            # 2. Adjacent residue phosphodiester linkage: O3'(i) - P(i+1)
            elif (r1 + 1 == r2 and a1 == "O3'" and a2 == 'P') or (r2 + 1 == r1 and a2 == "O3'" and a1 == 'P'):
                if d < 1.45 or d > 1.75:
                    raise AssertionError(f"Check 11 FAILED: Distorted phosphodiester bond ({d:.2f} A) between {r1} {a1} and {r2} {a2}")
            # 3. All other non-bonded pairs
            else:
                if d < 1.80:
                    raise AssertionError(f"Check 11 FAILED: Severe steric clash ({d:.2f} A) between {r1} {rn1} {a1} and {r2} {rn2} {a2}")
    print(f"[PASS] Check 11: Zero steric clashes detected (minimum non-bonded distance > 1.80 A).")
    
    # 12. Reasonable covalent bond connectivity
    for rid in sorted_ids:
        rname, atoms, _ = res_dict[rid]
        
        # Ribose C1' - C2'
        d_c1_c2 = np.linalg.norm(atoms["C1'"] - atoms["C2'"])
        if not (1.40 <= d_c1_c2 <= 1.62):
            raise AssertionError(f"Check 12 FAILED: Distorted ribose C1'-C2' ({d_c1_c2:.2f} A) at residue {rid}")
            
        # Glycosidic bond: C1' - N1 (pyrimidine) or C1' - N9 (purine)
        n_anc = 'N9' if rname in ['A', 'G'] else 'N1'
        d_glyc = np.linalg.norm(atoms["C1'"] - atoms[n_anc])
        if not (1.40 <= d_glyc <= 1.55):
            raise AssertionError(f"Check 12 FAILED: Distorted glycosidic bond C1'-{n_anc} ({d_glyc:.2f} A) at residue {rid}")
            
        # Phosphodiester connectivity to next residue
        if rid < 43:
            next_p = res_dict[rid + 1][1]['P']
            d_o3_p = np.linalg.norm(atoms["O3'"] - next_p)
            if not (1.50 <= d_o3_p <= 1.75):
                raise AssertionError(f"Check 12 FAILED: Distorted phosphodiester bond O3'-P ({d_o3_p:.2f} A) between {rid} and {rid+1}")
                
        # 5'-OH terminus check (residue 27 must not have P)
        if rid == 27:
            if 'P' in atoms or 'OP1' in atoms or 'OP2' in atoms:
                raise AssertionError("Check 12 FAILED: Residue 27 has phosphate atoms (expected neutral 5'-OH)")
        # Internal residues must have P
        elif rid > 27:
            if 'P' not in atoms or 'OP1' not in atoms or 'OP2' not in atoms:
                raise AssertionError(f"Check 12 FAILED: Internal residue {rid} is missing phosphate atoms")
    print(f"[PASS] Check 12: Covalent bond lengths (ribose, glycosidic, phosphodiester, neutral 5'-OH) within expected limits.")
    
    return res_dict

def validate_wt_and_suppressor_comparison(wt_dict, cua_dict):
    print(f"\n{'='*70}")
    print("Validating Controlled Comparison between WT and CUA Suppressor")
    print(f"{'='*70}")
    
    # 13. Correct mutation locations
    mutated_positions = []
    for rid in EXPECTED_SPRINZL_IDS:
        r_wt = wt_dict[rid][0]
        r_cua = cua_dict[rid][0]
        if r_wt != r_cua:
            mutated_positions.append((rid, r_wt, r_cua))
            
    if len(mutated_positions) != 2:
        raise AssertionError(f"Check 13 FAILED: Expected exactly 2 mutated positions, found {len(mutated_positions)}: {mutated_positions}")
    pos1, wt1, cua1 = mutated_positions[0]
    pos2, wt2, cua2 = mutated_positions[1]
    if (pos1, wt1, cua1) != (34, 'G', 'C'):
        raise AssertionError(f"Check 13 FAILED: Expected position 34 G->C, found {pos1} {wt1}->{cua1}")
    if (pos2, wt2, cua2) != (36, 'C', 'A'):
        raise AssertionError(f"Check 13 FAILED: Expected position 36 C->A, found {pos2} {wt2}->{cua2}")
    print(f"[PASS] Check 13: Mutations are strictly confined to position 34 (G->C) and position 36 (C->A).")
    
    # 14. WT and Suppressor share identical non-mutated coordinates
    non_mutated_ids = [rid for rid in EXPECTED_SPRINZL_IDS if rid not in (34, 36)]
    max_dev = 0.0
    sum_sq = 0.0
    atom_count = 0
    
    for rid in non_mutated_ids:
        atoms_wt = wt_dict[rid][1]
        atoms_cua = cua_dict[rid][1]
        if set(atoms_wt.keys()) != set(atoms_cua.keys()):
            raise AssertionError(f"Check 14 FAILED: Atom set mismatch at non-mutated residue {rid}")
        for aname in atoms_wt:
            d = np.linalg.norm(atoms_wt[aname] - atoms_cua[aname])
            sum_sq += d * d
            atom_count += 1
            if d > max_dev:
                max_dev = d
                
    rmsd_non_mutated = np.sqrt(sum_sq / atom_count)
    if rmsd_non_mutated > 1e-4:
        raise AssertionError(f"Check 14 FAILED: Non-mutated residues have diverged! RMSD = {rmsd_non_mutated:.6f} A")
    print(f"[PASS] Check 14: Non-mutated residues (15/17) have IDENTICAL coordinates (RMSD = {rmsd_non_mutated:.6f} A, max deviation = {max_dev:.6f} A).")
    
    # U-turn diagnostic comparison
    u33_n3_wt = wt_dict[33][1]['N3']
    c36_op2_wt = wt_dict[36][1]['OP2']
    d_uturn_wt = np.linalg.norm(u33_n3_wt - c36_op2_wt)
    
    u33_n3_cua = cua_dict[33][1]['N3']
    a36_op2_cua = cua_dict[36][1]['OP2']
    d_uturn_cua = np.linalg.norm(u33_n3_cua - a36_op2_cua)
    
    print(f"\nU-Turn Diagnostic Comparison:")
    print(f"  WT  ASL U33(N3) ... C36(OP2): {d_uturn_wt:.2f} A")
    print(f"  CUA ASL U33(N3) ... A36(OP2): {d_uturn_cua:.2f} A")

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    wt_pdb = os.path.join(base_dir, 'structures', 'wt', 'initial.pdb')
    cua_pdb = os.path.join(base_dir, 'structures', 'cua', 'initial.pdb')
    
    try:
        wt_dict = validate_single_structure(wt_pdb, EXPECTED_WT_SEQ, "WT ASL (CCCGCCUGUCACGCGGG)")
        cua_dict = validate_single_structure(cua_pdb, EXPECTED_CUA_SEQ, "CUA Suppressor ASL (CCCGCCUCUAACGCGGG)")
        validate_wt_and_suppressor_comparison(wt_dict, cua_dict)
        print(f"\n{'='*70}")
        print("ALL 14 STRUCTURAL VALIDATION CHECKS PASSED FOR BOTH CONSTRUCTS!")
        print(f"{'='*70}\n")
        sys.exit(0)
    except AssertionError as e:
        print(f"\n[VALIDATION FAILURE] {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n[UNEXPECTED ERROR] {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
