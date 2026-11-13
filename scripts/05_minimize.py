#!/usr/bin/env python3
"""
scripts/05_minimize.py - Phase 4B Energy Minimization Pipeline

Scientific Objective:
Perform energy minimization on the solvated Wild-Type and CUA Suppressor 17-nt
ASL explicit-solvent systems prepared in Phase 4A.

Purpose:
Remove unfavorable van der Waals steric clashes, relax strained solvent-solute
contacts, and optimize initial electrostatics prior to restrained equilibration.

Protocol Specifications (SDR-07, SDR-09, docs/05_simulation_protocol.md):
  - Force field: AMBER OL3 RNA (ff99bsc0chiOL3) + TIP3P + Joung-Cheatham KCl
  - Input: Serialized OpenMM system (simulations/{construct}/system/system.xml)
           and solvated coordinates (simulations/{construct}/system/solvated.pdb)
  - Two-Stage Minimization:
      Stage 1: Restrained Solute Minimization (RNA heavy atoms harmonically restrained
               with k = 5.0 kcal/(mol*A^2) ~ 2092 kJ/(mol*nm^2), max 2500 iterations)
      Stage 2: Full System Minimization (unrestrained, tolerance = 10.0 kJ/(mol*nm),
               max 5000 iterations)
  - Output:
      simulations/{construct}/minimized.pdb
      simulations/{construct}/minimization.log
      simulations/{construct}/minimization_metadata.yaml
      simulations/{construct}/minimization_metadata.json

STRICT PHASE BOUNDARY:
This script performs energy minimization ONLY.
No NVT/NPT equilibration, production MD, trajectory generation, RMSD/RMSF
analysis, or biological interpretation is performed in Phase 4B.
Energy minimization removes unfavorable local geometries but does not establish
conformational stability, folding, or biological function.
"""

import os
import sys
import time
import json
import yaml
import argparse
import hashlib
import datetime
import numpy as np

import openmm as mm
import openmm.app as app
import openmm.unit as unit

# Authoritative construct sequences and metadata
CONSTRUCT_INFO = {
    'wt': {
        'construct_id': 'wt_asl',
        'construct_name': 'Wild-Type ASL',
        'target_trna': 'tRNA-Asp-GTC-2-1',
        'sequence_5_to_3': 'CCCGCCUGUCACGCGGG',
        'expected_atoms': 15693,
        'expected_rna_atoms': 544,
        'expected_rna_heavy_atoms': 357,
        'expected_waters': 5035,
        'expected_k_ions': 30,
        'expected_cl_ions': 14,
        'dir_name': 'wt'
    },
    'cua': {
        'construct_id': 'suppressor_asl',
        'construct_name': 'CUA Suppressor ASL',
        'target_trna': 'tRNA-Asp-GTC-2-1 (CUA mutant)',
        'sequence_5_to_3': 'CCCGCCUCUAACGCGGG',
        'expected_atoms': 15686,
        'expected_rna_atoms': 543,
        'expected_rna_heavy_atoms': 356,
        'expected_waters': 5033,
        'expected_k_ions': 30,
        'expected_cl_ions': 14,
        'dir_name': 'cua'
    }
}

def compute_sha256(filepath):
    """Calculates SHA-256 checksum of a file."""
    h = hashlib.sha256()
    with open(filepath, 'rb') as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def select_best_platform(preferred='auto'):
    """Selects and configures the best available OpenMM platform."""
    available = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    
    if preferred != 'auto':
        if preferred not in available:
            raise ValueError(f"Requested platform '{preferred}' not available. Choose from {available}")
        plat = mm.Platform.getPlatformByName(preferred)
        props = {}
        if preferred in ['OpenCL', 'CUDA']:
            props['Precision'] = 'mixed'
        return plat, preferred, props
        
    # Auto-detection priority: CUDA > OpenCL > CPU > Reference
    for name in ['CUDA', 'OpenCL', 'CPU', 'Reference']:
        if name in available:
            plat = mm.Platform.getPlatformByName(name)
            props = {}
            if name in ['OpenCL', 'CUDA']:
                props['Precision'] = 'mixed'
            return plat, name, props
            
    raise RuntimeError("No compatible OpenMM platform discovered.")

def get_platform_hardware_info(platform, platform_name):
    """Queries device name and properties from an OpenMM platform."""
    info = {'platform_name': platform_name}
    try:
        # Create a minimal context to query properties
        s = mm.System()
        s.addParticle(1.0)
        integ = mm.VerletIntegrator(0.001)
        c = mm.Context(s, integ, platform)
        prop_names = platform.getPropertyNames()
        for p in prop_names:
            try:
                info[p] = platform.getPropertyValue(c, p)
            except Exception:
                pass
        del c
    except Exception as e:
        info['query_warning'] = str(e)
    return info

def extract_rna_heavy_atom_indices(topology):
    """Returns list of atom indices for all RNA heavy (non-hydrogen) atoms."""
    indices = []
    for res in topology.residues():
        if res.name in ['A', 'C', 'G', 'U']:
            for atom in res.atoms():
                if not atom.name.startswith('H'):
                    indices.append(atom.index)
    return indices

def validate_minimized_structure(pdb_path, initial_pdb_path, construct_info, initial_energy, final_energy):
    """
    Performs comprehensive post-minimization scientific and geometric validation.
    """
    pdb = app.PDBFile(pdb_path)
    init_pdb = app.PDBFile(initial_pdb_path)
    
    results = {}
    passed = True
    errors = []
    
    # Check 1: Total atom count
    n_atoms = len(pdb.positions)
    expected_atoms = construct_info['expected_atoms']
    if n_atoms != expected_atoms:
        passed = False
        errors.append(f"Atom count mismatch: found {n_atoms}, expected {expected_atoms}")
    results['total_atoms'] = n_atoms
    
    # Check 2: Coordinate validity (no NaN, Inf)
    coords = np.array([p.value_in_unit(unit.nanometer) for p in pdb.positions])
    if np.isnan(coords).any() or np.isinf(coords).any():
        passed = False
        errors.append("NaN or Inf coordinates detected in minimized structure!")
    results['nan_or_inf_detected'] = bool(np.isnan(coords).any() or np.isinf(coords).any())
    
    # Check 3: RNA residue count and sequence
    rna_residues = [r for r in pdb.topology.residues() if r.name in ['A', 'C', 'G', 'U']]
    seq = "".join(r.name for r in rna_residues)
    expected_seq = construct_info['sequence_5_to_3']
    if len(rna_residues) != 17:
        passed = False
        errors.append(f"Expected 17 RNA residues, found {len(rna_residues)}")
    if seq != expected_seq:
        passed = False
        errors.append(f"RNA sequence mismatch: found {seq}, expected {expected_seq}")
    results['rna_residues'] = len(rna_residues)
    results['rna_sequence'] = seq
    
    # Check 4: Residue numbering and Sprinzl mapping
    res_ids = [int(r.id) for r in rna_residues]
    if res_ids != list(range(1, 18)) and res_ids != list(range(27, 44)):
        passed = False
        errors.append(f"Residue numbering mismatch: found {res_ids}, expected [1..17] or [27..43]")
    results['residue_id_range'] = f"{res_ids[0]}-{res_ids[-1]}"
    results['sprinzl_mapping'] = "Residues 1-17 map 1-to-1 to Sprinzl positions 27-43"
    
    # Check 5: Terminal chemistry
    first_atoms = {a.name for a in rna_residues[0].atoms()}
    last_atoms = {a.name for a in rna_residues[-1].atoms()}
    has_ho5 = "HO5'" in first_atoms
    no_p = "P" not in first_atoms
    has_ho3 = "HO3'" in last_atoms
    if not (has_ho5 and no_p and has_ho3):
        passed = False
        errors.append("Terminal neutral 5'-OH / 3'-OH chemistry compromised during minimization!")
    results['terminal_chemistry_intact'] = bool(has_ho5 and no_p and has_ho3)
    
    # Check 6: Solvent and ion inventories
    waters = sum(1 for r in pdb.topology.residues() if r.name == 'HOH')
    k_ions = sum(1 for r in pdb.topology.residues() if r.name == 'K')
    cl_ions = sum(1 for r in pdb.topology.residues() if r.name == 'CL')
    if waters != construct_info['expected_waters']:
        passed = False
        errors.append(f"Water count mismatch: found {waters}, expected {construct_info['expected_waters']}")
    if k_ions != construct_info['expected_k_ions']:
        passed = False
        errors.append(f"K+ count mismatch: found {k_ions}, expected {construct_info['expected_k_ions']}")
    if cl_ions != construct_info['expected_cl_ions']:
        passed = False
        errors.append(f"Cl- count mismatch: found {cl_ions}, expected {construct_info['expected_cl_ions']}")
    results['water_molecules'] = waters
    results['k_ions'] = k_ions
    results['cl_ions'] = cl_ions
    
    # Check 7: Energy behavior (E_final <= E_initial)
    delta_e = final_energy - initial_energy
    if delta_e > 1e-3:
        passed = False
        errors.append(f"Minimization failed to reduce potential energy: delta_E = {delta_e:.3f} kJ/mol")
    results['initial_potential_energy_kj_per_mol'] = float(initial_energy)
    results['final_potential_energy_kj_per_mol'] = float(final_energy)
    results['energy_change_kj_per_mol'] = float(delta_e)
    
    # Check 8: RNA heavy-atom RMSD (structural preservation)
    rna_heavy_indices = extract_rna_heavy_atom_indices(pdb.topology)
    init_coords = np.array([p.value_in_unit(unit.nanometer) for p in init_pdb.positions])
    rna_init = init_coords[rna_heavy_indices]
    rna_final = coords[rna_heavy_indices]
    rna_rmsd_angstrom = float(np.sqrt(np.mean(np.sum((rna_final - rna_init)**2, axis=1))) * 10.0)
    results['rna_heavy_atom_rmsd_angstrom'] = rna_rmsd_angstrom
    if rna_rmsd_angstrom > 2.0:
        passed = False
        errors.append(f"Excessive RNA deformation during minimization: RMSD = {rna_rmsd_angstrom:.3f} A > 2.0 A")
        
    # Check 9: Non-bonded steric clash check (subsample RNA to find minimum distance)
    # Minimum interatomic distance within RNA heavy atoms
    diffs = rna_final[:, np.newaxis, :] - rna_final[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diffs**2, axis=-1)) * 10.0 # to Angstroms
    np.fill_diagonal(dists, 999.0)
    min_dist_angstrom = float(np.min(dists))
    results['minimum_rna_heavy_distance_angstrom'] = min_dist_angstrom
    if min_dist_angstrom < 1.0:
        passed = False
        errors.append(f"Severe steric clash detected after minimization: minimum distance = {min_dist_angstrom:.3f} A < 1.0 A")
        
    results['passed'] = passed
    results['errors'] = errors
    return results

def run_minimization(construct_key, config, platform_pref='auto', no_restrained=False,
                     override_tol=None, override_max_iter=None, base_dir=None):
    """
    Executes the complete energy minimization workflow for a single construct.
    """
    if base_dir is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
    info = CONSTRUCT_INFO[construct_key]
    construct_dir = os.path.join(base_dir, 'simulations', info['dir_name'])
    system_dir = os.path.join(construct_dir, 'system')
    
    xml_path = os.path.join(system_dir, 'system.xml')
    solvated_pdb_path = os.path.join(system_dir, 'solvated.pdb')
    
    if not os.path.exists(xml_path):
        raise FileNotFoundError(f"Missing Phase 4A system XML: {xml_path}")
    if not os.path.exists(solvated_pdb_path):
        raise FileNotFoundError(f"Missing Phase 4A solvated PDB: {solvated_pdb_path}")
        
    output_pdb_path = os.path.join(construct_dir, 'minimized.pdb')
    log_path = os.path.join(construct_dir, 'minimization.log')
    metadata_yaml_path = os.path.join(construct_dir, 'minimization_metadata.yaml')
    metadata_json_path = os.path.join(construct_dir, 'minimization_metadata.json')
    
    # Read minimization parameters from config
    min_cfg = config.get('minimization', {})
    tol_kj_mol_nm = override_tol if override_tol is not None else min_cfg.get('tolerance_kj_per_mol_nm', 10.0)
    max_iter_unrestrained = override_max_iter if override_max_iter is not None else min_cfg.get('max_iterations_unrestrained', 5000)
    
    stage1_cfg = min_cfg.get('stage1_restrained', {})
    stage1_enabled = False if no_restrained else stage1_cfg.get('enabled', True)
    k_kcal_a2 = stage1_cfg.get('restraint_k_kcal_per_mol_a2', 5.0)
    max_iter_restrained = stage1_cfg.get('max_iterations_restrained', 2500)
    
    start_time = datetime.datetime.now(datetime.timezone.utc)
    
    with open(log_path, 'w', encoding='utf-8') as log_file:
        def log(msg):
            t_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            line = f"[{t_str}] {msg}"
            print(line)
            log_file.write(line + "\n")
            log_file.flush()
            
        log("=" * 80)
        log(f"Phase 4B Energy Minimization: {info['construct_name']} ({info['target_trna']})")
        log("=" * 80)
        log(f"Input System XML: {xml_path}")
        log(f"Input Solvated PDB: {solvated_pdb_path}")
        log(f"Output Minimized PDB: {output_pdb_path}")
        log(f"Protocol Parameters:")
        log(f"  - Stage 1 (Restrained Solute): {'ENABLED' if stage1_enabled else 'DISABLED'}")
        if stage1_enabled:
            log(f"    * Restraint constant (k): {k_kcal_a2} kcal/(mol*A^2)")
            log(f"    * Max iterations: {max_iter_restrained}")
        log(f"  - Stage 2 (Full System): ENABLED")
        log(f"    * Convergence tolerance: {tol_kj_mol_nm} kJ/(mol*nm)")
        log(f"    * Max iterations: {max_iter_unrestrained}")
        
        # 1. Load System XML
        log("Deserializing Phase 4A OpenMM System XML (no reparameterization)...")
        with open(xml_path, 'r', encoding='utf-8') as f:
            system = mm.XmlSerializer.deserialize(f.read())
        log(f"System loaded successfully: {system.getNumParticles()} particles, {system.getNumForces()} force groups.")
        
        # 2. Load Solvated PDB
        log("Loading Phase 4A Solvated PDB coordinates...")
        pdb = app.PDBFile(solvated_pdb_path)
        log(f"PDB loaded successfully: {len(pdb.positions)} coordinates.")
        
        if system.getNumParticles() != len(pdb.positions):
            raise AssertionError(f"Particle count mismatch: XML ({system.getNumParticles()}) vs PDB ({len(pdb.positions)})")
            
        # 3. Setup Platform
        log(f"Selecting OpenMM computation platform (preference: {platform_pref})...")
        platform, plat_name, plat_props = select_best_platform(platform_pref)
        for k_prop, v_prop in plat_props.items():
            platform.setPropertyDefaultValue(k_prop, v_prop)
        hw_info = get_platform_hardware_info(platform, plat_name)
        log(f"Configured platform: {plat_name} with properties {plat_props}")
        if 'DeviceName' in hw_info:
            log(f"Hardware compute device: {hw_info['DeviceName']}")
            
        # 4. Harmonic Position Restraints (Stage 1)
        rna_heavy_indices = extract_rna_heavy_atom_indices(pdb.topology)
        log(f"Identified {len(rna_heavy_indices)} RNA heavy atoms for harmonic position restraints.")
        
        k_kj_nm2 = (k_kcal_a2 * unit.kilocalories_per_mole / unit.angstrom**2).value_in_unit(
            unit.kilojoules_per_mole / unit.nanometer**2
        )
        
        # Add CustomExternalForce for restraints
        restraint_force = mm.CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
        restraint_force.addGlobalParameter("k", 0.0) # Start with k=0 for initial energy evaluation
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")
        
        for idx in rna_heavy_indices:
            pos_nm = pdb.positions[idx].value_in_unit(unit.nanometer)
            restraint_force.addParticle(idx, [pos_nm[0], pos_nm[1], pos_nm[2]])
            
        system.addForce(restraint_force)
        
        # 5. Build Simulation Context
        integrator = mm.LangevinMiddleIntegrator(310.0*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtosecond)
        sim = app.Simulation(pdb.topology, system, integrator, platform)
        sim.context.setPositions(pdb.positions)
        
        # Verify box vectors
        box_vectors = pdb.topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            sim.context.setPeriodicBoxVectors(*box_vectors)
            
        # 6. Evaluate Initial Unrestrained Potential Energy
        sim.context.setParameter("k", 0.0)
        state_init = sim.context.getState(getEnergy=True)
        e_initial = state_init.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        log(f"Initial Phase 4A Solvated Potential Energy (unrestrained): {e_initial:14.3f} kJ/mol")
        
        # 7. Execute Stage 1: Restrained Solute Minimization
        e_stage1 = e_initial
        if stage1_enabled:
            log("Executing Stage 1: Restrained Solute Minimization (solvent/ion relaxation)...")
            log(f"  Applying RNA heavy-atom harmonic restraints: k = {k_kj_nm2:.2f} kJ/(mol*nm^2)")
            sim.context.setParameter("k", k_kj_nm2)
            
            t_s1_start = time.time()
            sim.minimizeEnergy(
                tolerance=tol_kj_mol_nm * unit.kilojoules_per_mole / unit.nanometer,
                maxIterations=max_iter_restrained
            )
            t_s1_elapsed = time.time() - t_s1_start
            
            # Record post-stage 1 unrestrained energy
            sim.context.setParameter("k", 0.0)
            state_s1 = sim.context.getState(getEnergy=True)
            e_stage1 = state_s1.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            delta_s1 = e_stage1 - e_initial
            log(f"Stage 1 completed in {t_s1_elapsed:.2f}s:")
            log(f"  Post-Stage 1 Potential Energy (unrestrained): {e_stage1:14.3f} kJ/mol")
            log(f"  Stage 1 Energy Drop (delta_E):                {delta_s1:14.3f} kJ/mol")
        else:
            log("Stage 1 skipped by configuration (pure single-stage unrestrained minimization).")
            
        # 8. Execute Stage 2: Full System Minimization
        log("Executing Stage 2: Full System Minimization (unrestrained)...")
        sim.context.setParameter("k", 0.0) # Ensure zero restraint
        
        t_s2_start = time.time()
        sim.minimizeEnergy(
            tolerance=tol_kj_mol_nm * unit.kilojoules_per_mole / unit.nanometer,
            maxIterations=max_iter_unrestrained
        )
        t_s2_elapsed = time.time() - t_s2_start
        
        state_final = sim.context.getState(getEnergy=True, getPositions=True)
        e_final = state_final.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        delta_total = e_final - e_initial
        delta_s2 = e_final - e_stage1
        
        log(f"Stage 2 completed in {t_s2_elapsed:.2f}s:")
        log(f"  Final System Potential Energy:                 {e_final:14.3f} kJ/mol")
        log(f"  Stage 2 Additional Drop:                      {delta_s2:14.3f} kJ/mol")
        log(f"  TOTAL MINIMIZATION ENERGY DROP (delta_E):      {delta_total:14.3f} kJ/mol")
        
        # 9. Save Minimized PDB
        log(f"Writing minimized coordinates to: {output_pdb_path}...")
        with open(output_pdb_path, 'w', encoding='utf-8') as f_pdb:
            app.PDBFile.writeFile(pdb.topology, state_final.getPositions(), f_pdb, keepIds=True)
        pdb_sha = compute_sha256(output_pdb_path)
        log(f"Minimized PDB written successfully (SHA-256: {pdb_sha})")
        
        # 10. Automated Validation
        log("Running post-minimization automated validation...")
        val_results = validate_minimized_structure(
            output_pdb_path, solvated_pdb_path, info, e_initial, e_final
        )
        for k_res, v_res in val_results.items():
            if k_res not in ['passed', 'errors']:
                log(f"  * {k_res}: {v_res}")
        if val_results['passed']:
            log("[PASS] All post-minimization structural, chemical, and energetic checks verified!")
        else:
            log(f"[FAIL] Post-minimization validation failed: {val_results['errors']}")
            raise AssertionError(f"Validation failed for {info['construct_name']}: {val_results['errors']}")
            
        end_time = datetime.datetime.now(datetime.timezone.utc)
        total_walltime = (end_time - start_time).total_seconds()
        
        # 11. Compile Metadata
        metadata = {
            'construct': construct_key,
            'construct_id': info['construct_id'],
            'construct_name': info['construct_name'],
            'selected_trna_identifier': info['target_trna'],
            'asl_sequence_5_to_3': info['sequence_5_to_3'],
            'input_files': {
                'solvated_pdb': solvated_pdb_path,
                'solvated_pdb_sha256': compute_sha256(solvated_pdb_path),
                'system_xml': xml_path,
                'system_xml_sha256': compute_sha256(xml_path)
            },
            'output_files': {
                'minimized_pdb': output_pdb_path,
                'minimized_pdb_sha256': pdb_sha,
                'minimization_log': log_path
            },
            'openmm_environment': {
                'openmm_version': mm.__version__,
                'platform_used': plat_name,
                'platform_hardware': hw_info,
                'precision_mode': plat_props.get('Precision', 'default'),
                'python_version': sys.version.split()[0]
            },
            'minimization_protocol': {
                'tolerance_kj_per_mol_nm': float(tol_kj_mol_nm),
                'tolerance_kcal_per_mol_a': float(tol_kj_mol_nm * 0.02390057),
                'max_iterations_unrestrained': int(max_iter_unrestrained),
                'stage1_restrained': {
                    'enabled': bool(stage1_enabled),
                    'restraint_k_kcal_per_mol_a2': float(k_kcal_a2) if stage1_enabled else None,
                    'restraint_k_kj_per_mol_nm2': float(k_kj_nm2) if stage1_enabled else None,
                    'max_iterations': int(max_iter_restrained) if stage1_enabled else None,
                    'restrained_atom_count': len(rna_heavy_indices) if stage1_enabled else 0
                }
            },
            'energy_results': {
                'initial_potential_energy_kj_per_mol': float(e_initial),
                'stage1_unrestrained_energy_kj_per_mol': float(e_stage1),
                'final_potential_energy_kj_per_mol': float(e_final),
                'delta_energy_kj_per_mol': float(delta_total),
                'energy_units': 'kJ/mol',
                'energy_decreased': bool(delta_total <= 0.0)
            },
            'validation_checks': val_results,
            'timing': {
                'start_time_utc': start_time.isoformat(),
                'end_time_utc': end_time.isoformat(),
                'walltime_seconds': float(total_walltime)
            },
            'status': 'PASS',
            'disclaimer': (
                "Energy minimization removes unfavorable local steric contacts and relaxes solvent-solute "
                "geometry, but does NOT establish conformational stability, native folding, or biological function."
            )
        }
        
        # Save Metadata YAML and JSON
        with open(metadata_yaml_path, 'w', encoding='utf-8') as f_yaml:
            yaml.dump(metadata, f_yaml, sort_keys=False, default_flow_style=False)
        with open(metadata_json_path, 'w', encoding='utf-8') as f_json:
            json.dump(metadata, f_json, indent=2)
            
        log(f"Metadata exported to: {metadata_yaml_path}")
        log(f"Metadata exported to: {metadata_json_path}")
        log(f"Phase 4B Minimization finished successfully for {info['construct_name']} (Verdict: PASS).\n")
        
    return metadata

def main():
    parser = argparse.ArgumentParser(
        description="Phase 4B: Energy Minimization of WT and CUA Solvated ASL Systems"
    )
    parser.add_argument(
        '--construct', choices=['wt', 'cua', 'all'], default='all',
        help="Construct to minimize ('wt', 'cua', or 'all' [default])"
    )
    parser.add_argument(
        '--config', default='config/simulation_config.yaml',
        help="Path to simulation config YAML (default: config/simulation_config.yaml)"
    )
    parser.add_argument(
        '--platform', default='auto',
        choices=['auto', 'OpenCL', 'CPU', 'CUDA', 'Reference'],
        help="OpenMM compute platform (default: auto)"
    )
    parser.add_argument(
        '--no-restrained-stage', action='store_true',
        help="Skip Stage 1 restrained solute minimization"
    )
    parser.add_argument(
        '--tolerance', type=float, default=None,
        help="Override convergence tolerance in kJ/(mol*nm)"
    )
    parser.add_argument(
        '--max-iterations', type=int, default=None,
        help="Override maximum iterations for unrestrained stage"
    )
    parser.add_argument(
        '--validate-only', action='store_true',
        help="Run post-minimization validation only on existing minimized.pdb files"
    )
    
    args = parser.parse_args()
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_path = os.path.join(base_dir, args.config)
    
    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
        
    constructs = ['wt', 'cua'] if args.construct == 'all' else [args.construct]
    
    if args.validate_only:
        print("=" * 80)
        print("PHASE 4B: RUNNING POST-MINIMIZATION VALIDATION ONLY")
        print("=" * 80)
        all_ok = True
        for c in constructs:
            info = CONSTRUCT_INFO[c]
            p_min = os.path.join(base_dir, 'simulations', info['dir_name'], 'minimized.pdb')
            p_solv = os.path.join(base_dir, 'simulations', info['dir_name'], 'system', 'solvated.pdb')
            p_meta = os.path.join(base_dir, 'simulations', info['dir_name'], 'minimization_metadata.json')
            
            if not os.path.exists(p_min) or not os.path.exists(p_meta):
                print(f"[ERROR] Missing minimized files for {c}. Run minimization first.")
                all_ok = False
                continue
            with open(p_meta, 'r') as f:
                meta = json.load(f)
            e_init = meta['energy_results']['initial_potential_energy_kj_per_mol']
            e_final = meta['energy_results']['final_potential_energy_kj_per_mol']
            v = validate_minimized_structure(p_min, p_solv, info, e_init, e_final)
            print(f"\nValidation for {info['construct_name']}: {'PASS' if v['passed'] else 'FAIL'}")
            for k, val in v.items():
                if k not in ['passed', 'errors']:
                    print(f"  * {k}: {val}")
            if not v['passed']:
                print(f"  ERRORS: {v['errors']}")
                all_ok = False
        sys.exit(0 if all_ok else 1)
        
    print("=" * 80)
    print("PHASE 4B: SOLVATED ASL ENERGY MINIMIZATION")
    print(f"Constructs: {constructs} | Platform preference: {args.platform}")
    print("=" * 80)
    
    results = {}
    for c in constructs:
        results[c] = run_minimization(
            construct_key=c,
            config=config,
            platform_pref=args.platform,
            no_restrained=args.no_restrained_stage,
            override_tol=args.tolerance,
            override_max_iter=args.max_iterations,
            base_dir=base_dir
        )
        
    # Comparative summary
    if len(constructs) == 2:
        wt_m = results['wt']
        cua_m = results['cua']
        print("\n" + "=" * 80)
        print("PHASE 4B COMPARATIVE MINIMIZATION SUMMARY (WT vs. CUA)")
        print("=" * 80)
        print(f"{'Metric':<34} | {'Wild-Type (WT)':<20} | {'Suppressor (CUA)':<20}")
        print("-" * 80)
        print(f"{'Sequence (5\'->3\')':<34} | {wt_m['asl_sequence_5_to_3']:<20} | {cua_m['asl_sequence_5_to_3']:<20}")
        print(f"{'Total System Particles':<34} | {wt_m['validation_checks']['total_atoms']:<20} | {cua_m['validation_checks']['total_atoms']:<20}")
        print(f"{'Platform / Precision':<34} | {wt_m['openmm_environment']['platform_used']}/{wt_m['openmm_environment']['precision_mode']:<14} | {cua_m['openmm_environment']['platform_used']}/{cua_m['openmm_environment']['precision_mode']:<14}")
        print(f"{'Initial Energy (kJ/mol)':<34} | {wt_m['energy_results']['initial_potential_energy_kj_per_mol']:<20.3f} | {cua_m['energy_results']['initial_potential_energy_kj_per_mol']:<20.3f}")
        print(f"{'Post-Stage 1 Energy (kJ/mol)':<34} | {wt_m['energy_results']['stage1_unrestrained_energy_kj_per_mol']:<20.3f} | {cua_m['energy_results']['stage1_unrestrained_energy_kj_per_mol']:<20.3f}")
        print(f"{'Final Energy (kJ/mol)':<34} | {wt_m['energy_results']['final_potential_energy_kj_per_mol']:<20.3f} | {cua_m['energy_results']['final_potential_energy_kj_per_mol']:<20.3f}")
        print(f"{'Total Energy Drop delta_E (kJ/mol)':<34} | {wt_m['energy_results']['delta_energy_kj_per_mol']:<20.3f} | {cua_m['energy_results']['delta_energy_kj_per_mol']:<20.3f}")
        print(f"{'RNA Heavy-Atom RMSD (A)':<34} | {wt_m['validation_checks']['rna_heavy_atom_rmsd_angstrom']:<20.4f} | {cua_m['validation_checks']['rna_heavy_atom_rmsd_angstrom']:<20.4f}")
        print(f"{'Min RNA Distance (A)':<34} | {wt_m['validation_checks']['minimum_rna_heavy_distance_angstrom']:<20.4f} | {cua_m['validation_checks']['minimum_rna_heavy_distance_angstrom']:<20.4f}")
        print(f"{'Minimization Status':<34} | {wt_m['status']:<20} | {cua_m['status']:<20}")
        print("=" * 80)
        print("\n[PHASE 4B COMPLETE] Energy minimization finished and verified for all constructs.")
        print("IMPORTANT: No equilibration or MD has been performed. Ready for Phase 4C.\n")

if __name__ == '__main__':
    main()
