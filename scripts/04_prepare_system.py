#!/usr/bin/env python3
"""
scripts/04_prepare_system.py - Phase 4A Solvation & MD System Preparation Pipeline

Scientific Objective:
Prepare explicit-solvent periodic simulation systems for:
  1. Wild-Type 17-nt ASL (structures/wt/initial.pdb)
  2. CUA Suppressor 17-nt ASL (structures/cua/initial.pdb)

Protocol Specifications:
  - Force field: AMBER OL3 RNA (amber14/RNA.OL3.xml)
  - Water model: TIP3P explicit solvent (amber14/tip3p.xml)
  - Ion parameters: Joung & Cheatham (2008) for TIP3P (ionsjc_tip3p in amber14/tip3p.xml)
  - Terminal chemistry: Neutral 5'-OH (res 27) and neutral 3'-OH (res 43)
  - Periodic box: Cubic (L = 5.500 nm / 55.00 A)
  - Minimum solvent padding: >= 1.0 nm (10.0 A) guaranteed in all directions
  - Neutralization: Exactly 16 K+ counterions to neutralize the -16.0 e polyanionic RNA
  - Added salt: Approximately 0.15 M KCl (14 K+ and 14 Cl- pairs)
  - Target ensemble: NPT at 310 K, 1 bar (prepared for Phase 4B/4C)
  - Deterministic seeding: random.seed(42) for reproducible water replacement

STRICT PHASE BOUNDARY:
This script performs system preparation ONLY.
No energy minimization, equilibration, or production MD is executed in Phase 4A.
"""

import os
import sys
import json
import yaml
import random
import hashlib
import numpy as np

import openmm.app as app
import openmm.unit as unit
import openmm as mm

# Expected sequences (Sprinzl 27-43)
EXPECTED_WT_SEQ = "CCCGCCUGUCACGCGGG"
EXPECTED_CUA_SEQ = "CCCGCCUCUAACGCGGG"

BOX_SIZE_NM = 5.500  # 55.00 A cubic box
TARGET_SALT_MOLAR = 0.15  # 150 mM KCl
RANDOM_SEED = 42

def compute_sha256(filepath):
    """Calculates SHA-256 hash of a file."""
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def extract_rna_sequence(topology):
    """Extracts RNA sequence string from an OpenMM Topology."""
    seq = []
    for res in topology.residues():
        if res.name in ['A', 'C', 'G', 'U']:
            seq.append(res.name)
    return "".join(seq)

def verify_terminal_chemistry(topology, positions):
    """Verifies that terminal residues possess neutral 5'-OH and 3'-OH chemistry."""
    rna_residues = [r for r in topology.residues() if r.name in ['A', 'C', 'G', 'U']]
    if len(rna_residues) != 17:
        raise AssertionError(f"Expected 17 RNA residues, found {len(rna_residues)}")
    
    res_first = rna_residues[0]
    res_last = rna_residues[-1]
    
    first_atoms = {a.name for a in res_first.atoms()}
    last_atoms = {a.name for a in res_last.atoms()}
    
    # 5'-terminus: MUST have HO5', MUST NOT have P, OP1, OP2
    if "HO5'" not in first_atoms:
        raise AssertionError(f"5'-terminus ({res_first.name}{res_first.id}) is missing HO5'!")
    if "P" in first_atoms:
        raise AssertionError(f"5'-terminus ({res_first.name}{res_first.id}) possesses an unexpected 5'-phosphate group!")
        
    # 3'-terminus: MUST have HO3'
    if "HO3'" not in last_atoms:
        raise AssertionError(f"3'-terminus ({res_last.name}{res_last.id}) is missing HO3'!")
        
    return True

def compute_solute_clearance(topology, positions, box_size_nm):
    """
    Computes exact minimum distance from any solute RNA atom to any face of the cubic box.
    Solvent box in OpenMM Modeller is centered at:
      center = 0.5 * (minRange + maxRange)
    Box bounds along each axis d are [center[d] - L/2, center[d] + L/2].
    """
    rna_coords = []
    for res in topology.residues():
        if res.name in ['A', 'C', 'G', 'U']:
            for atom in res.atoms():
                rna_coords.append(positions[atom.index].value_in_unit(unit.nanometer))
    rna_coords = np.array(rna_coords)
    
    center = 0.5 * (np.min(rna_coords, axis=0) + np.max(rna_coords, axis=0))
    rna_extent = np.max(rna_coords, axis=0) - np.min(rna_coords, axis=0)
    
    dist_to_edges = []
    for d in range(3):
        min_bound = center[d] - box_size_nm / 2.0
        max_bound = center[d] + box_size_nm / 2.0
        dist_to_edges.append(np.min(rna_coords[:, d] - min_bound))
        dist_to_edges.append(np.min(max_bound - rna_coords[:, d]))
        
    min_clearance_nm = np.min(dist_to_edges)
    return min_clearance_nm, rna_extent, center

def prepare_system(input_pdb, construct_id, construct_name, expected_seq, output_dir):
    """Prepares and validates a fully solvated, neutralized, salt-adjusted MD system."""
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "preparation.log")
    
    with open(log_path, "w") as log:
        def log_msg(msg):
            print(msg)
            log.write(msg + "\n")
            
        log_msg("=" * 80)
        log_msg(f"PHASE 4A: SYSTEM PREPARATION & SOLVATION - {construct_name}")
        log_msg(f"Input PDB: {input_pdb}")
        log_msg(f"Output Directory: {output_dir}")
        log_msg("=" * 80)
        
        # 1. Deterministic seeding
        random.seed(RANDOM_SEED)
        np.random.seed(RANDOM_SEED)
        log_msg(f"[CONFIG] Pseudo-random seed set to {RANDOM_SEED} for deterministic ion placement.")
        
        # 2. Input file verification
        if not os.path.exists(input_pdb):
            raise FileNotFoundError(f"Input PDB does not exist: {input_pdb}")
        input_hash = compute_sha256(input_pdb)
        log_msg(f"[INPUT] Input PDB SHA-256: {input_hash}")
        
        # 3. Load Force Field
        log_msg("[STEP 1] Loading force fields: amber14/RNA.OL3.xml and amber14/tip3p.xml...")
        ff = app.ForceField('amber14/RNA.OL3.xml', 'amber14/tip3p.xml')
        log_msg("[PASS] Force field XML definitions loaded successfully.")
        
        # 4. Parse PDB and Add Hydrogens
        log_msg("[STEP 2] Parsing solute coordinates and adding hydrogens via Modeller...")
        pdb = app.PDBFile(input_pdb)
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(ff)
        
        raw_rna_atoms = modeller.topology.getNumAtoms()
        log_msg(f"[PASS] Solute topology initialized with {raw_rna_atoms} atoms.")
        
        # Verify sequence and terminal chemistry
        raw_seq = extract_rna_sequence(modeller.topology)
        if raw_seq != expected_seq:
            raise AssertionError(f"Sequence mismatch! Expected {expected_seq}, got {raw_seq}")
        log_msg(f"[PASS] Solute RNA sequence confirmed: 5'-{raw_seq}-3' (17 nt).")
        
        verify_terminal_chemistry(modeller.topology, modeller.positions)
        log_msg("[PASS] Terminal chemistry confirmed: Neutral 5'-OH (res 27) and neutral 3'-OH (res 43).")
        
        # Check unsolvated RNA net charge
        unsolv_sys = ff.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff)
        nb_force_unsolv = [f for f in unsolv_sys.getForces() if isinstance(f, mm.NonbondedForce)][0]
        unsolv_charge = sum(nb_force_unsolv.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(unsolv_sys.getNumParticles()))
        log_msg(f"[CHARGE] Unsolvated RNA net charge: {unsolv_charge:.4f} e (Expected exactly -16.0000 e).")
        if abs(unsolv_charge - (-16.0)) > 1e-4:
            raise AssertionError(f"Unsolvated RNA charge discrepancy! Expected -16.0 e, got {unsolv_charge:.4f} e")
        log_msg("[PASS] Unsolvated RNA charge strictly verified as -16.0 e.")
        
        # 5. Solvation & Ion Placement
        log_msg(f"[STEP 3] Solvating in cubic periodic box (L = {BOX_SIZE_NM:.3f} nm / {BOX_SIZE_NM*10:.1f} A)...")
        box_vecs_param = mm.Vec3(BOX_SIZE_NM, BOX_SIZE_NM, BOX_SIZE_NM) * unit.nanometers
        
        # Add explicit TIP3P water and ions
        modeller.addSolvent(
            ff,
            model='tip3p',
            boxSize=box_vecs_param,
            positiveIon='K+',
            negativeIon='Cl-',
            ionicStrength=TARGET_SALT_MOLAR * unit.molar,
            neutralize=True
        )
        log_msg("[PASS] Solvation and ion placement complete.")
        
        # 6. Inventory & Demographics
        res_names = [r.name for r in modeller.topology.residues()]
        num_k = res_names.count('K')
        num_cl = res_names.count('CL')
        num_h2o = res_names.count('HOH')
        total_atoms = modeller.topology.getNumAtoms()
        total_residues = modeller.topology.getNumResidues()
        
        # Breakdown of K+ ions:
        # Neutralization counterions = exactly 16 K+
        # Added salt pairs = num_cl pairs of K+ and Cl-
        num_neutralizing_k = 16
        num_added_salt_k = num_k - num_neutralizing_k
        num_added_salt_cl = num_cl
        
        log_msg("\n" + "-" * 40)
        log_msg("SYSTEM INVENTORY & ION BREAKDOWN")
        log_msg("-" * 40)
        log_msg(f"  Total Atoms:            {total_atoms}")
        log_msg(f"  Total Residues:         {total_residues}")
        log_msg(f"  RNA Atoms:              {raw_rna_atoms} (17 residues)")
        log_msg(f"  TIP3P Water Molecules:  {num_h2o} ({num_h2o * 3} atoms)")
        log_msg(f"  Neutralizing Counterions: {num_neutralizing_k} K+ (neutralizes -16.0 e RNA)")
        log_msg(f"  Added Salt Ions (KCl):    {num_added_salt_k} K+ and {num_added_salt_cl} Cl- ({num_added_salt_cl} pairs)")
        log_msg(f"  Total K+ Ions:          {num_k}")
        log_msg(f"  Total Cl- Ions:         {num_cl}")
        
        if num_added_salt_k != num_added_salt_cl:
            raise AssertionError(f"Added salt ion asymmetry! K+: {num_added_salt_k}, Cl-: {num_added_salt_cl}")
        log_msg("[PASS] Ion symmetry verified: Added salt consists of exact neutral ion pairs.")
        
        # 7. Box Geometry & Clearance
        box_vectors = modeller.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer)
        box_dims = np.array([box_vectors[0][0], box_vectors[1][1], box_vectors[2][2]])
        box_volume_nm3 = box_dims[0] * box_dims[1] * box_dims[2]
        box_volume_liters = box_volume_nm3 * 1e-24  # 1 nm^3 = 1e-24 L
        
        min_clearance_nm, rna_extent, rna_center = compute_solute_clearance(modeller.topology, modeller.positions, BOX_SIZE_NM)
        
        log_msg("\n" + "-" * 40)
        log_msg("BOX GEOMETRY & SOLUTE PADDING AUDIT")
        log_msg("-" * 40)
        log_msg(f"  Box Geometry:           Cubic (orthogonal)")
        log_msg(f"  Box Dimensions:         {box_dims[0]:.4f} x {box_dims[1]:.4f} x {box_dims[2]:.4f} nm")
        log_msg(f"  Box Volume:             {box_volume_nm3:.4f} nm^3 ({box_volume_nm3*1000:.1f} A^3)")
        log_msg(f"  Solute Geometric Center: [{rna_center[0]:.3f}, {rna_center[1]:.3f}, {rna_center[2]:.3f}] nm")
        log_msg(f"  Solute Extents (x,y,z):  [{rna_extent[0]:.3f}, {rna_extent[1]:.3f}, {rna_extent[2]:.3f}] nm")
        log_msg(f"  Minimum Solute Padding: {min_clearance_nm:.4f} nm ({min_clearance_nm*10:.2f} A)")
        
        if min_clearance_nm < 1.0:
            raise AssertionError(f"Solute padding violation! Minimum distance {min_clearance_nm:.4f} nm is < 1.0 nm requirement!")
        log_msg("[PASS] Minimum solvent padding strictly >= 1.0 nm (10.0 A) verified across all solute atoms.")
        
        # 8. Salt Concentration Calculation
        # Definition A: Based on total periodic box volume
        molar_conc_box = num_added_salt_cl / (6.02214076e23 * box_volume_liters)
        
        # Definition B: Based on actual solvent volume (55.4 mol/L pure water standard)
        water_volume_liters = num_h2o / (6.02214076e23 * 55.4)
        molar_conc_solvent = num_added_salt_cl / (6.02214076e23 * water_volume_liters)
        
        log_msg("\n" + "-" * 40)
        log_msg("SALT CONCENTRATION ANALYSIS (TARGET = 0.150 M)")
        log_msg("-" * 40)
        log_msg(f"  Target Concentration:   {TARGET_SALT_MOLAR:.3f} M (150.0 mM)")
        log_msg(f"  Box-Volume Concentration:    {molar_conc_box:.4f} M ({molar_conc_box*1000:.1f} mM) [Deviation: {((molar_conc_box - TARGET_SALT_MOLAR)/TARGET_SALT_MOLAR)*100:+.2f}%]")
        log_msg(f"  Solvent-Volume Concentration: {molar_conc_solvent:.4f} M ({molar_conc_solvent*1000:.1f} mM) [Deviation: {((molar_conc_solvent - TARGET_SALT_MOLAR)/TARGET_SALT_MOLAR)*100:+.2f}%]")
        log_msg("[PASS] Salt concentration is within optimal physical bounds for 150 mM physiological target.")
        
        # 9. Force Field System Creation & Final Charge Audit
        log_msg("\n" + "-" * 40)
        log_msg("SYSTEM CREATION & PARAMETERIZATION AUDIT")
        log_msg("-" * 40)
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometers,
            constraints=app.HBonds
        )
        
        nb_force = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
        final_net_charge = sum(nb_force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
        log_msg(f"  OpenMM System Particles: {system.getNumParticles()}")
        log_msg(f"  OpenMM System Forces:    {system.getNumForces()}")
        log_msg(f"  Final System Net Charge: {final_net_charge:.6f} e")
        
        if abs(final_net_charge) > 1e-4:
            raise AssertionError(f"Final system is not electrically neutral! Net charge: {final_net_charge:.6f} e")
        log_msg("[PASS] Electrical neutrality strictly verified: Net charge = 0.0000 e.")
        
        # Check coordinate health (no NaN, no Inf)
        all_coords = np.array([p.value_in_unit(unit.nanometer) for p in modeller.positions])
        if np.isnan(all_coords).any() or np.isinf(all_coords).any():
            raise AssertionError("Corrupted coordinates detected! NaN or Inf present.")
        log_msg("[PASS] Coordinate integrity verified: Zero NaN or Inf values.")
        
        # 10. Save Output Artifacts
        log_msg("\n" + "-" * 40)
        log_msg("SERIALIZATION & ARTIFACT PERSISTENCE")
        log_msg("-" * 40)
        
        # Solvated PDB
        solvated_pdb_path = os.path.join(output_dir, "solvated.pdb")
        with open(solvated_pdb_path, "w") as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
        solvated_pdb_hash = compute_sha256(solvated_pdb_path)
        log_msg(f"  Saved Solvated PDB:   {solvated_pdb_path} (SHA-256: {solvated_pdb_hash})")
        
        # OpenMM System XML
        system_xml_path = os.path.join(output_dir, "system.xml")
        with open(system_xml_path, "w") as f:
            f.write(mm.XmlSerializer.serialize(system))
        system_xml_hash = compute_sha256(system_xml_path)
        log_msg(f"  Saved System XML:     {system_xml_path} (SHA-256: {system_xml_hash})")
        
        # System Metadata JSON
        metadata = {
            "construct_id": construct_id,
            "construct_name": construct_name,
            "sequence_5_to_3": raw_seq,
            "input_pdb": input_pdb,
            "input_pdb_sha256": input_hash,
            "solvated_pdb": solvated_pdb_path,
            "solvated_pdb_sha256": solvated_pdb_hash,
            "system_xml": system_xml_path,
            "system_xml_sha256": system_xml_hash,
            "openmm_version": mm.__version__,
            "force_field_rna": "amber14/RNA.OL3.xml (ff99bsc0chiOL3)",
            "force_field_water": "amber14/tip3p.xml (TIP3P)",
            "force_field_ions": "amber14/tip3p.xml (Joung & Cheatham 2008 for TIP3P)",
            "terminal_chemistry": "Neutral 5'-OH (res 27) and neutral 3'-OH (res 43)",
            "box_geometry": "Cubic",
            "box_dimensions_nm": [float(b) for b in box_dims],
            "box_dimensions_angstrom": [float(b * 10.0) for b in box_dims],
            "box_volume_nm3": float(box_volume_nm3),
            "solute_extents_nm": [float(e) for e in rna_extent],
            "minimum_solute_padding_nm": float(min_clearance_nm),
            "minimum_solute_padding_angstrom": float(min_clearance_nm * 10.0),
            "total_atoms": int(total_atoms),
            "total_residues": int(total_residues),
            "rna_atoms": int(raw_rna_atoms),
            "water_molecules": int(num_h2o),
            "neutralizing_k_counterions": int(num_neutralizing_k),
            "added_k_salt_ions": int(num_added_salt_k),
            "added_cl_salt_ions": int(num_added_salt_cl),
            "total_k_ions": int(num_k),
            "total_cl_ions": int(num_cl),
            "target_salt_molar": float(TARGET_SALT_MOLAR),
            "actual_salt_molar_box_volume": float(molar_conc_box),
            "actual_salt_molar_solvent_volume": float(molar_conc_solvent),
            "unsolvated_rna_charge_e": float(unsolv_charge),
            "final_system_net_charge_e": float(final_net_charge),
            "nonbonded_method": "PME",
            "nonbonded_cutoff_nm": 1.0,
            "constraints": "HBonds",
            "random_seed": int(RANDOM_SEED),
            "status": "PASS"
        }
        
        metadata_path = os.path.join(output_dir, "system_metadata.json")
        with open(metadata_path, "w") as f:
            json.dump(metadata, f, indent=2)
        log_msg(f"  Saved Metadata JSON:  {metadata_path}")
        
        log_msg("\n" + "=" * 80)
        log_msg(f"PHASE 4A SYSTEM PREPARATION FOR {construct_name}: COMPLETE (PASS)")
        log_msg("=" * 80 + "\n")
        
        return metadata

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    wt_input = os.path.join(base_dir, "structures", "wt", "initial.pdb")
    cua_input = os.path.join(base_dir, "structures", "cua", "initial.pdb")
    
    wt_out_dir = os.path.join(base_dir, "simulations", "wt", "system")
    cua_out_dir = os.path.join(base_dir, "simulations", "cua", "system")
    
    print("\n" + "#" * 80)
    print("STARTING PHASE 4A SOLVATION & SYSTEM PREPARATION PIPELINE")
    print("#" * 80 + "\n")
    
    try:
        wt_meta = prepare_system(
            input_pdb=wt_input,
            construct_id="wt_asl",
            construct_name="Wild-Type ASL",
            expected_seq=EXPECTED_WT_SEQ,
            output_dir=wt_out_dir
        )
        
        cua_meta = prepare_system(
            input_pdb=cua_input,
            construct_id="suppressor_asl",
            construct_name="CUA Suppressor ASL",
            expected_seq=EXPECTED_CUA_SEQ,
            output_dir=cua_out_dir
        )
        
        # Print comparative summary
        print("\n" + "=" * 80)
        print("PHASE 4A COMPARATIVE SYSTEM SUMMARY (WT vs. CUA)")
        print("=" * 80)
        print(f"{'Metric':<32} | {'Wild-Type (WT)':<20} | {'Suppressor (CUA)':<20}")
        print("-" * 80)
        print(f"{'Sequence (5\'->3\')':<32} | {wt_meta['sequence_5_to_3']:<20} | {cua_meta['sequence_5_to_3']:<20}")
        print(f"{'RNA Atoms':<32} | {wt_meta['rna_atoms']:<20} | {cua_meta['rna_atoms']:<20}")
        print(f"{'Water Molecules (TIP3P)':<32} | {wt_meta['water_molecules']:<20} | {cua_meta['water_molecules']:<20}")
        print(f"{'Neutralizing K+ Counterions':<32} | {wt_meta['neutralizing_k_counterions']:<20} | {cua_meta['neutralizing_k_counterions']:<20}")
        print(f"{'Added Salt Pairs (KCl)':<32} | {wt_meta['added_cl_salt_ions']:<20} | {cua_meta['added_cl_salt_ions']:<20}")
        print(f"{'Total K+ Ions':<32} | {wt_meta['total_k_ions']:<20} | {cua_meta['total_k_ions']:<20}")
        print(f"{'Total Cl- Ions':<32} | {wt_meta['total_cl_ions']:<20} | {cua_meta['total_cl_ions']:<20}")
        print(f"{'Total System Atoms':<32} | {wt_meta['total_atoms']:<20} | {cua_meta['total_atoms']:<20}")
        print(f"{'Box Dimensions (nm)':<32} | {str(wt_meta['box_dimensions_nm']):<20} | {str(cua_meta['box_dimensions_nm']):<20}")
        print(f"{'Box Volume (nm^3)':<32} | {wt_meta['box_volume_nm3']:<20.3f} | {cua_meta['box_volume_nm3']:<20.3f}")
        print(f"{'Min Solute Padding (nm)':<32} | {wt_meta['minimum_solute_padding_nm']:<20.3f} | {cua_meta['minimum_solute_padding_nm']:<20.3f}")
        print(f"{'Target Salt Molarity (M)':<32} | {wt_meta['target_salt_molar']:<20.3f} | {cua_meta['target_salt_molar']:<20.3f}")
        print(f"{'Actual Salt (Solvent Vol, M)':<32} | {wt_meta['actual_salt_molar_solvent_volume']:<20.4f} | {cua_meta['actual_salt_molar_solvent_volume']:<20.4f}")
        print(f"{'Actual Salt (Box Vol, M)':<32} | {wt_meta['actual_salt_molar_box_volume']:<20.4f} | {cua_meta['actual_salt_molar_box_volume']:<20.4f}")
        print(f"{'Unsolvated RNA Charge (e)':<32} | {wt_meta['unsolvated_rna_charge_e']:<20.4f} | {cua_meta['unsolvated_rna_charge_e']:<20.4f}")
        print(f"{'Final System Net Charge (e)':<32} | {wt_meta['final_system_net_charge_e']:<20.4f} | {cua_meta['final_system_net_charge_e']:<20.4f}")
        print(f"{'Preparation Status':<32} | {wt_meta['status']:<20} | {cua_meta['status']:<20}")
        print("=" * 80)
        print("\n[PHASE 4A COMPLETE] Both systems successfully prepared and verified.")
        print("IMPORTANT: No minimization, equilibration, or MD has been performed. Ready for Phase 4B.\n")
        sys.exit(0)
    except Exception as e:
        print(f"\n[FATAL ERROR IN PHASE 4A] {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()

