"""
Phase 6: Trajectory Preparation and PBC Imaging.
Extracts solute-only (17-nt ASL) analysis-ready trajectories with explicit PBC imaging.
"""
import os
import sys
import json
import mdtraj as md

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

REPLICATES = [
    ('wt', 1), ('wt', 2), ('wt', 3),
    ('cua', 1), ('cua', 2), ('cua', 3)
]

def extract_solute_trajectories():
    print("=" * 80)
    print("PHASE 6: PREPARING ANALYSIS-READY PBC-IMAGED SOLUTE TRAJECTORIES")
    print("=" * 80)
    
    metadata = {}
    
    for construct, rep_num in REPLICATES:
        rep_tag = f"{construct}_rep{rep_num}"
        rep_dir = os.path.join(BASE_DIR, 'simulations', rep_tag)
        dcd_path = os.path.join(rep_dir, 'trajectory', 'trajectory.dcd')
        pdb_path = os.path.join(rep_dir, 'final_state', 'final_structure.pdb')
        
        solute_dcd = os.path.join(rep_dir, 'trajectory', 'solute_imaged.dcd')
        solute_pdb = os.path.join(rep_dir, 'trajectory', 'solute_ref.pdb')
        
        print(f"\nProcessing {rep_tag}...")
        print(f"  Source DCD: {dcd_path}")
        print(f"  Source PDB: {pdb_path}")
        
        # Load trajectory
        traj = md.load(dcd_path, top=pdb_path)
        print(f"  Loaded {traj.n_frames} frames with {traj.n_atoms} atoms.")
        
        # PBC Treatment: image whole molecules
        print("  Applying PBC imaging (image_molecules)...")
        traj_imaged = traj.image_molecules()
        
        # Select solute atoms
        solute_idx = traj.topology.select('not water and not (resname K or resname CL)')
        print(f"  Selected {len(solute_idx)} solute RNA atoms (17 nucleotides).")
        
        # Atom slice
        solute_traj = traj_imaged.atom_slice(solute_idx)
        
        # Save analysis-ready reference PDB and DCD
        solute_traj[0].save_pdb(solute_pdb)
        solute_traj.save_dcd(solute_dcd)
        sz_mb = os.path.getsize(solute_dcd) / (1024 * 1024)
        print(f"  Saved analysis-ready DCD: {solute_dcd} ({sz_mb:.2f} MB)")
        print(f"  Saved reference solute PDB: {solute_pdb}")
        
        metadata[rep_tag] = {
            'construct': construct,
            'replicate': rep_num,
            'source_dcd': os.path.relpath(dcd_path, BASE_DIR),
            'source_pdb': os.path.relpath(pdb_path, BASE_DIR),
            'solute_dcd': os.path.relpath(solute_dcd, BASE_DIR),
            'solute_pdb': os.path.relpath(solute_pdb, BASE_DIR),
            'frames': solute_traj.n_frames,
            'solute_atoms': solute_traj.n_atoms,
            'pbc_treatment': 'mdtraj.image_molecules()',
            'stride_ps': 10.0,
            'duration_ns': 20.0
        }
        
    meta_path = os.path.join(BASE_DIR, 'analysis', 'trajectory_preparation_metadata.json')
    with open(meta_path, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2)
    print(f"\nMetadata written to: {meta_path}")
    print("Trajectory extraction and PBC imaging complete for all 6 replicates.\n")

if __name__ == '__main__':
    extract_solute_trajectories()
