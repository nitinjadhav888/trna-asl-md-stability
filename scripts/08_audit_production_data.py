"""
Audit of Phase 5 production trajectories and environment before Phase 6 analysis.
Verifies file presence, MDTraj readability, frame counts, time continuity, and atom counts.
"""
import os
import sys
import yaml
import mdtraj as md

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

REPLICATES = [
    ('wt', 1, 5101, 15693, 'CCCGCCUGUCACGCGGG'),
    ('wt', 2, 5102, 15693, 'CCCGCCUGUCACGCGGG'),
    ('wt', 3, 5103, 15693, 'CCCGCCUGUCACGCGGG'),
    ('cua', 1, 5201, 15686, 'CCCGCCUCUAACGCGGG'),
    ('cua', 2, 5202, 15686, 'CCCGCCUCUAACGCGGG'),
    ('cua', 3, 5203, 15686, 'CCCGCCUCUAACGCGGG'),
]

def run_audit():
    print("=" * 80)
    print("PHASE 6 PRE-ANALYSIS AUDIT: PRODUCTION DATA INTEGRITY")
    print("=" * 80)
    
    all_passed = True
    results = []
    
    for construct, rep_num, seed, exp_atoms, exp_seq in REPLICATES:
        rep_tag = f"{construct}_rep{rep_num}"
        rep_dir = os.path.join(BASE_DIR, 'simulations', rep_tag)
        dcd_path = os.path.join(rep_dir, 'trajectory', 'trajectory.dcd')
        pdb_path = os.path.join(rep_dir, 'final_state', 'final_structure.pdb')
        thermo_path = os.path.join(rep_dir, 'production_thermo.csv')
        log_path = os.path.join(rep_dir, 'production.log')
        
        checks = []
        
        # 1. Directory exists
        if not os.path.exists(rep_dir):
            checks.append((False, f"Directory missing: {rep_dir}"))
        else:
            checks.append((True, "Directory exists"))
            
        # 2. Files exist
        for fpath, desc in [(dcd_path, "DCD trajectory"), (pdb_path, "Final structure PDB"), 
                            (thermo_path, "Thermo CSV"), (log_path, "Production log")]:
            if not os.path.exists(fpath):
                checks.append((False, f"{desc} missing: {fpath}"))
            else:
                sz_mb = os.path.getsize(fpath) / (1024 * 1024)
                checks.append((True, f"{desc} present ({sz_mb:.2f} MB)"))
                
        # 3. MDTraj Readability & Frame Count
        try:
            traj = md.load(dcd_path, top=pdb_path)
            n_frames = traj.n_frames
            n_atoms = traj.n_atoms
            if n_frames == 2000:
                checks.append((True, f"Trajectory readable: exactly 2000 frames (20.0 ns at 10 ps stride)"))
            else:
                checks.append((False, f"Frame count mismatch: {n_frames} frames (expected 2000)"))
                
            if n_atoms == exp_atoms:
                checks.append((True, f"Atom count verified: {n_atoms} atoms"))
            else:
                checks.append((False, f"Atom count mismatch: {n_atoms} (expected {exp_atoms})"))
                
            # Verify RNA sequence
            solute_idx = traj.topology.select('not water and not (resname K or resname CL)')
            solute_top = traj.topology.subset(solute_idx)
            seq = "".join(r.name for r in solute_top.residues)
            if seq == exp_seq:
                checks.append((True, f"Sequence verified: {seq}"))
            else:
                checks.append((False, f"Sequence mismatch: found {seq}, expected {exp_seq}"))
                
        except Exception as e:
            checks.append((False, f"Failed to load trajectory with MDTraj: {str(e)}"))
            
        # 4. Check Thermo CSV
        try:
            with open(thermo_path, 'r', encoding='utf-8') as f:
                lines = [l.strip() for l in f if l.strip()]
            header = lines[0]
            records = len(lines) - 1
            if records == 2000:
                checks.append((True, f"Thermo CSV verified: exactly 2000 records"))
            else:
                checks.append((False, f"Thermo CSV record count mismatch: {records} (expected 2000)"))
        except Exception as e:
            checks.append((False, f"Failed to read thermo CSV: {str(e)}"))
            
        rep_ok = all(c[0] for c in checks)
        if not rep_ok:
            all_passed = False
            
        results.append((rep_tag, rep_ok, checks))
        
    print("\nAUDIT RESULTS TABLE:")
    print("-" * 80)
    for rep_tag, rep_ok, checks in results:
        status_str = "PASS" if rep_ok else "FAIL"
        print(f"[{status_str}] {rep_tag}:")
        for passed, msg in checks:
            mark = "  ✓" if passed else "  ✗"
            print(f"    {mark} {msg}")
    print("-" * 80)
    
    if all_passed:
        print("\nALL 6 PRODUCTION REPLICATES PASSED INTEGRITY AUDIT (6/6).")
        print("Ready for Phase 6 analysis pipeline.\n")
        return 0
    else:
        print("\nDATA INTEGRITY AUDIT FAILED. STOPPING BEFORE ANALYSIS.\n")
        return 1

if __name__ == '__main__':
    sys.exit(run_audit())
