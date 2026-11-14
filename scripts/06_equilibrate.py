#!/usr/bin/env python3
"""
scripts/06_equilibrate.py - Phase 4C Restrained NVT/NPT Equilibration Pipeline

Scientific Objective:
Perform staged NVT heating and NPT density equilibration with progressive
harmonic restraint release on the minimized Wild-Type and CUA Suppressor 17-nt
ASL explicit-solvent systems from Phase 4B.

Thermodynamic Target Conditions:
  - Temperature: 310.0 K (physiological)
  - Pressure: 1.0 bar (1.01325 bar, isotropic Monte Carlo barostat)
  - Ensemble: NVT heating (200 ps) followed by NPT density equilibration (1000 ps)
  - Force field: AMBER OL3 RNA + TIP3P water + Joung-Cheatham KCl
  - Restraint schedule (on RNA heavy atoms):
      NVT:            k = 5.0 kcal/(mol*A^2) (ramp 10 -> 310 K over 150 ps, 50 ps plateau)
      NPT Stage 1:    k = 2.0 kcal/(mol*A^2) (250 ps)
      NPT Stage 2:    k = 1.0 kcal/(mol*A^2) (250 ps)
      NPT Stage 3:    k = 0.2 kcal/(mol*A^2) (250 ps)
      NPT Stage 4:    k = 0.0 kcal/(mol*A^2) (250 ps, unrestrained pre-production relaxation)
  - Total equilibration time: 1.2 ns (600,000 steps at 2.0 fs)

STRICT PHASE BOUNDARY:
This script performs equilibration ONLY.
No production MD, replicate sampling, trajectory analysis, RMSD/RMSF
analysis, or biological interpretation is performed in Phase 4C.
Equilibration prepares the molecular system for production dynamics by relaxing
solvent density and heating the system to target thermodynamic conditions,
but does NOT establish conformational stability, native folding, or biological function.
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

# Unit conversion constants
KCAL_TO_KJ = 4.184
KCAL_PER_MOL_A2_TO_KJ_PER_MOL_NM2 = 418.4
BOLTZMANN_CONSTANT_R = 8.314462618e-3  # kJ / (mol * K)
DA_PER_NM3_TO_G_PER_CM3 = 1.66053906660e-3

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
        s = mm.System()
        s.addParticle(1.0)
        integ = mm.VerletIntegrator(0.001)
        c = mm.Context(s, integ, platform)
        for p in platform.getPropertyNames():
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

def compute_instantaneous_temperature(ke_kj_mol, system):
    """Calculates instantaneous temperature from kinetic energy and degrees of freedom."""
    dof = 3 * system.getNumParticles() - system.getNumConstraints() - 3
    if dof <= 0:
        return 0.0
    return (2.0 * ke_kj_mol) / (dof * BOLTZMANN_CONSTANT_R)

def compute_box_density(system, box_vectors):
    """Calculates box volume (nm^3) and physical density (g/cm^3)."""
    total_mass_da = sum(system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(system.getNumParticles()))
    vol_nm3 = box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]
    density_g_cm3 = (total_mass_da / vol_nm3) * DA_PER_NM3_TO_G_PER_CM3
    return float(vol_nm3), float(density_g_cm3)

def validate_equilibrated_system(pdb_path, minimized_pdb_path, construct_info, npt_thermo_csv):
    """
    Performs automated scientific and physical validation of the equilibrated system.
    """
    pdb = app.PDBFile(pdb_path)
    min_pdb = app.PDBFile(minimized_pdb_path)
    
    results = {}
    passed = True
    errors = []
    
    # 1. Total atom count
    n_atoms = len(pdb.positions)
    expected_atoms = construct_info['expected_atoms']
    if n_atoms != expected_atoms:
        passed = False
        errors.append(f"Atom count mismatch: found {n_atoms}, expected {expected_atoms}")
    results['total_atoms'] = n_atoms
    
    # 2. Coordinate validity (no NaN/Inf)
    coords = np.array([p.value_in_unit(unit.nanometer) for p in pdb.positions])
    if np.isnan(coords).any() or np.isinf(coords).any():
        passed = False
        errors.append("NaN or Inf detected in final coordinates!")
    results['nan_or_inf_detected'] = bool(np.isnan(coords).any() or np.isinf(coords).any())
    
    # 3. RNA residue count and sequence
    rna_residues = [r for r in pdb.topology.residues() if r.name in ['A', 'C', 'G', 'U']]
    seq = "".join(r.name for r in rna_residues)
    expected_seq = construct_info['sequence_5_to_3']
    if len(rna_residues) != 17:
        passed = False
        errors.append(f"Expected 17 RNA residues, found {len(rna_residues)}")
    if seq != expected_seq:
        passed = False
        errors.append(f"Sequence mismatch: expected {expected_seq}, found {seq}")
    results['rna_residues'] = len(rna_residues)
    results['rna_sequence'] = seq
    
    # 4. Terminal chemistry
    first_atoms = {a.name for a in rna_residues[0].atoms()}
    last_atoms = {a.name for a in rna_residues[-1].atoms()}
    has_ho5 = "HO5'" in first_atoms
    no_p = "P" not in first_atoms
    has_ho3 = "HO3'" in last_atoms
    if not (has_ho5 and no_p and has_ho3):
        passed = False
        errors.append("Terminal 5'-OH / 3'-OH chemistry disrupted during equilibration!")
    results['terminal_chemistry_intact'] = bool(has_ho5 and no_p and has_ho3)
    
    # 5. Solvent and ion inventories
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
    
    # 6. RNA heavy-atom RMSD relative to minimized structure (Kabsch superposition)
    rna_heavy_indices = extract_rna_heavy_atom_indices(pdb.topology)
    min_coords = np.array([p.value_in_unit(unit.angstrom) for p in min_pdb.positions])
    eq_coords = np.array([p.value_in_unit(unit.angstrom) for p in pdb.positions])
    rna_min = min_coords[rna_heavy_indices]
    rna_eq = eq_coords[rna_heavy_indices]
    
    raw_rmsd = float(np.sqrt(np.mean(np.sum((rna_eq - rna_min)**2, axis=1))))
    results['raw_displacement_rmsd_angstrom'] = raw_rmsd
    
    # Kabsch rotational and translational superposition
    c_min = rna_min - np.mean(rna_min, axis=0)
    c_eq = rna_eq - np.mean(rna_eq, axis=0)
    H = c_min.T @ c_eq
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    V = Vt.T
    if d < 0:
        V[:, -1] *= -1
    R = V @ U.T
    c_min_rot = c_min @ R.T
    
    kabsch_rmsd = float(np.sqrt(np.mean(np.sum((c_eq - c_min_rot)**2, axis=1))))
    results['rna_heavy_atom_rmsd_angstrom'] = kabsch_rmsd
    if kabsch_rmsd > 3.5:
        passed = False
        errors.append(f"Excessive RNA structural deformation during equilibration: Kabsch RMSD = {kabsch_rmsd:.3f} A > 3.5 A")
        
    # 7. Minimum interatomic distance in RNA
    diffs = rna_eq[:, np.newaxis, :] - rna_eq[np.newaxis, :, :]
    dists = np.sqrt(np.sum(diffs**2, axis=-1)) * 10.0
    np.fill_diagonal(dists, 999.0)
    min_dist = float(np.min(dists))
    results['minimum_rna_distance_angstrom'] = min_dist
    if min_dist < 1.0:
        passed = False
        errors.append(f"Steric clash in equilibrated RNA: min distance = {min_dist:.3f} A < 1.0 A")
        
    # 8. Thermodynamic statistics from NPT log
    if os.path.exists(npt_thermo_csv):
        import csv
        temps, densities, volumes, pes = [], [], [], []
        with open(npt_thermo_csv, 'r') as f_csv:
            reader = csv.DictReader(f_csv)
            for row in reader:
                temps.append(float(row['temp_k']))
                densities.append(float(row['density_g_cm3']))
                volumes.append(float(row['volume_nm3']))
                pes.append(float(row['pe_kj_mol']))
                
        # Analyze second half of NPT equilibration (last 50% of frames)
        half_idx = len(temps) // 2
        last_temps = temps[half_idx:]
        last_densities = densities[half_idx:]
        last_volumes = volumes[half_idx:]
        
        mean_t = float(np.mean(last_temps))
        std_t = float(np.std(last_temps))
        mean_density = float(np.mean(last_densities))
        std_density = float(np.std(last_densities))
        mean_vol = float(np.mean(last_volumes))
        
        results['npt_late_mean_temperature_k'] = mean_t
        results['npt_late_temperature_std_k'] = std_t
        results['npt_late_mean_density_g_cm3'] = mean_density
        results['npt_late_density_std_g_cm3'] = std_density
        results['npt_late_mean_volume_nm3'] = mean_vol
        
        # Check temperature within 310 +/- 5 K
        if abs(mean_t - 310.0) > 5.0:
            passed = False
            errors.append(f"Equilibrated temperature out of bounds: mean = {mean_t:.2f} K (expected 310 +/- 5 K)")
            
        # Check density within physical range for water/ions (0.95 to 1.05 g/cm3)
        if not (0.95 <= mean_density <= 1.05):
            passed = False
            errors.append(f"Equilibrated density unphysical: mean = {mean_density:.4f} g/cm^3")
            
    results['passed'] = passed
    results['errors'] = errors
    return results

def run_equilibration(construct_key, config, platform_pref='auto', scale_factor=1.0, base_dir=None):
    """
    Executes the full staged NVT/NPT equilibration workflow for a single construct.
    """
    if base_dir is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
    info = CONSTRUCT_INFO[construct_key]
    construct_dir = os.path.join(base_dir, 'simulations', info['dir_name'])
    eq_dir = os.path.join(construct_dir, 'equilibration')
    nvt_dir = os.path.join(eq_dir, 'nvt')
    npt_dir = os.path.join(eq_dir, 'npt')
    
    os.makedirs(nvt_dir, exist_ok=True)
    os.makedirs(npt_dir, exist_ok=True)
    
    # Input paths
    xml_path = os.path.join(construct_dir, 'system', 'system.xml')
    min_pdb_path = os.path.join(construct_dir, 'minimized.pdb')
    
    if not os.path.exists(xml_path):
        raise FileNotFoundError(f"Missing Phase 4A system XML: {xml_path}")
    if not os.path.exists(min_pdb_path):
        raise FileNotFoundError(f"Missing Phase 4B minimized PDB: {min_pdb_path}")
        
    # Output paths
    log_path = os.path.join(eq_dir, 'equilibration.log')
    metadata_yaml_path = os.path.join(eq_dir, 'equilibration_metadata.yaml')
    metadata_json_path = os.path.join(eq_dir, 'equilibration_metadata.json')
    final_pdb_path = os.path.join(eq_dir, 'final_equilibrated.pdb')
    final_state_path = os.path.join(eq_dir, 'final_state.xml')
    
    nvt_pdb_path = os.path.join(nvt_dir, 'nvt_equilibrated.pdb')
    nvt_state_path = os.path.join(nvt_dir, 'nvt_state.xml')
    nvt_csv_path = os.path.join(nvt_dir, 'nvt_thermo.csv')
    npt_csv_path = os.path.join(npt_dir, 'npt_thermo.csv')
    
    # Read equilibration parameters from config
    eq_cfg = config.get('equilibration', {})
    nvt_cfg = eq_cfg.get('nvt', {})
    npt_cfg = eq_cfg.get('npt', {})
    
    timestep_fs = eq_cfg.get('timestep_fs', 2.0)
    friction_ps = eq_cfg.get('friction_coeff_ps', 1.0)
    
    # Scale durations if requested
    nvt_dur_ps = nvt_cfg.get('duration_ps', 200.0) * scale_factor
    nvt_ramp_ps = nvt_cfg.get('ramp_duration_ps', 150.0) * scale_factor
    nvt_plateau_ps = nvt_cfg.get('plateau_duration_ps', 50.0) * scale_factor
    nvt_k_kcal = nvt_cfg.get('restraint_k_kcal_per_mol_a2', 5.0)
    
    npt_dur_ps = npt_cfg.get('duration_ps', 1000.0) * scale_factor
    npt_temp_k = npt_cfg.get('temperature_kelvin', 310.0)
    npt_press_bar = npt_cfg.get('pressure_bar', 1.0)
    barostat_freq = npt_cfg.get('barostat_frequency_steps', 25)
    
    npt_stages = npt_cfg.get('stages', [
        {'name': 'stage1', 'duration_ps': 250.0, 'restraint_k_kcal_per_mol_a2': 2.0},
        {'name': 'stage2', 'duration_ps': 250.0, 'restraint_k_kcal_per_mol_a2': 1.0},
        {'name': 'stage3', 'duration_ps': 250.0, 'restraint_k_kcal_per_mol_a2': 0.2},
        {'name': 'stage4', 'duration_ps': 250.0, 'restraint_k_kcal_per_mol_a2': 0.0},
    ])
    
    # Deterministic seed
    seeds = eq_cfg.get('random_seeds', {'wt': 4101, 'cua': 4201})
    random_seed = seeds.get(construct_key, 4101 if construct_key == 'wt' else 4201)
    
    report_interval_ps = 1.0 * (scale_factor if scale_factor < 1.0 else 1.0)
    report_steps = max(1, int(round((report_interval_ps * 1000.0) / timestep_fs)))
    
    start_time_utc = datetime.datetime.now(datetime.timezone.utc)
    
    with open(log_path, 'w', encoding='utf-8') as log_file:
        def log(msg):
            t_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            line = f"[{t_str}] {msg}"
            print(line)
            log_file.write(line + "\n")
            log_file.flush()
            
        log("=" * 80)
        log(f"Phase 4C Restrained Equilibration: {info['construct_name']} ({info['target_trna']})")
        log("=" * 80)
        log(f"Input System XML:       {xml_path}")
        log(f"Input Minimized PDB:    {min_pdb_path}")
        log(f"Equilibration Output:   {eq_dir}")
        log(f"Deterministic Seed:     {random_seed}")
        log(f"Scale Factor:           {scale_factor:.4f}")
        log(f"NVT Heating Duration:   {nvt_dur_ps:.1f} ps ({int(nvt_dur_ps*1000/timestep_fs)} steps)")
        log(f"NPT Density Duration:   {npt_dur_ps:.1f} ps ({int(npt_dur_ps*1000/timestep_fs)} steps)")
        log(f"Thermostat / Barostat:  LangevinMiddle (gamma={friction_ps}/ps) / MonteCarloBarostat ({npt_press_bar} bar, freq={barostat_freq})")
        
        # 1. Load System XML and Minimized PDB
        log("Deserializing Phase 4A OpenMM System XML...")
        with open(xml_path, 'r', encoding='utf-8') as f:
            system_nvt = mm.XmlSerializer.deserialize(f.read())
        pdb = app.PDBFile(min_pdb_path)
        log(f"Loaded system ({system_nvt.getNumParticles()} particles) and coordinates ({len(pdb.positions)} atoms).")
        
        # 2. Add Harmonic Restraints
        rna_heavy_indices = extract_rna_heavy_atom_indices(pdb.topology)
        log(f"Adding harmonic position restraints on {len(rna_heavy_indices)} RNA heavy atoms...")
        
        k_nvt_kj_nm2 = nvt_k_kcal * KCAL_PER_MOL_A2_TO_KJ_PER_MOL_NM2
        restraint_force = mm.CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
        restraint_force.addGlobalParameter("k", k_nvt_kj_nm2)
        restraint_force.addPerParticleParameter("x0")
        restraint_force.addPerParticleParameter("y0")
        restraint_force.addPerParticleParameter("z0")
        
        # Reference coordinates are the Phase 4B minimized positions
        for idx in rna_heavy_indices:
            pos_nm = pdb.positions[idx].value_in_unit(unit.nanometer)
            restraint_force.addParticle(idx, [pos_nm[0], pos_nm[1], pos_nm[2]])
            
        system_nvt.addForce(restraint_force)
        
        # 3. Setup Platform
        log(f"Selecting computation platform (preference: {platform_pref})...")
        platform, plat_name, plat_props = select_best_platform(platform_pref)
        for kp, vp in plat_props.items():
            platform.setPropertyDefaultValue(kp, vp)
        hw_info = get_platform_hardware_info(platform, plat_name)
        log(f"Platform selected: {plat_name} with properties {plat_props}")
        if 'DeviceName' in hw_info:
            log(f"Compute device: {hw_info['DeviceName']}")
            
        # ----------------------------------------------------------------------
        # PART I: NVT RESTRAINED HEATING (10 K -> 310 K)
        # ----------------------------------------------------------------------
        log("-" * 80)
        log("PART I: NVT Restrained Heating (10 K -> 310 K)")
        log("-" * 80)
        
        integrator_nvt = mm.LangevinMiddleIntegrator(
            10.0 * unit.kelvin,
            friction_ps / unit.picosecond,
            timestep_fs * unit.femtosecond
        )
        integrator_nvt.setRandomNumberSeed(random_seed)
        
        sim_nvt = app.Simulation(pdb.topology, system_nvt, integrator_nvt, platform)
        sim_nvt.context.setPositions(pdb.positions)
        
        box_vectors = pdb.topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            sim_nvt.context.setPeriodicBoxVectors(*box_vectors)
            
        log("Assigning initial velocities from Maxwell-Boltzmann distribution at 10.0 K...")
        sim_nvt.context.setVelocitiesToTemperature(10.0 * unit.kelvin, random_seed)
        
        total_nvt_steps = int(round((nvt_dur_ps * 1000.0) / timestep_fs))
        ramp_steps = int(round((nvt_ramp_ps * 1000.0) / timestep_fs))
        plateau_steps = total_nvt_steps - ramp_steps
        
        log(f"NVT Schedule: Ramp over {ramp_steps} steps ({nvt_ramp_ps:.1f} ps), then plateau for {plateau_steps} steps ({nvt_plateau_ps:.1f} ps).")
        
        # Prepare CSV logging for NVT
        nvt_csv_file = open(nvt_csv_path, 'w', encoding='utf-8')
        nvt_csv_file.write("step,time_ps,target_temp_k,temp_k,pe_kj_mol,ke_kj_mol,total_energy_kj_mol\n")
        
        step_chunks = report_steps
        current_nvt_step = 0
        t_nvt_start = time.time()
        
        while current_nvt_step < total_nvt_steps:
            steps_to_run = min(step_chunks, total_nvt_steps - current_nvt_step)
            
            # Update target temperature linearly during ramp
            if current_nvt_step < ramp_steps:
                frac = current_nvt_step / float(ramp_steps)
                target_T = 10.0 + (300.0 * frac)
            else:
                target_T = 310.0
                
            integrator_nvt.setTemperature(target_T * unit.kelvin)
            sim_nvt.step(steps_to_run)
            current_nvt_step += steps_to_run
            
            state = sim_nvt.context.getState(getEnergy=True)
            ke = state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            inst_T = compute_instantaneous_temperature(ke, system_nvt)
            time_ps = (current_nvt_step * timestep_fs) / 1000.0
            
            nvt_csv_file.write(f"{current_nvt_step},{time_ps:.3f},{target_T:.2f},{inst_T:.2f},{pe:.3f},{ke:.3f},{pe+ke:.3f}\n")
            
            # Print periodic log updates every ~20 ps or at end
            if current_nvt_step % (report_steps * 20) == 0 or current_nvt_step == total_nvt_steps:
                log(f"  NVT Step {current_nvt_step:6d}/{total_nvt_steps:6d} ({time_ps:5.1f} ps) | Target T: {target_T:5.1f} K | Inst T: {inst_T:5.1f} K | PE: {pe:11.1f} kJ/mol")
                
        nvt_csv_file.close()
        t_nvt_elapsed = time.time() - t_nvt_start
        log(f"NVT heating stage completed in {t_nvt_elapsed:.2f}s.")
        
        # Save NVT final state and coordinates
        nvt_state = sim_nvt.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
        with open(nvt_pdb_path, 'w', encoding='utf-8') as f_pdb:
            app.PDBFile.writeFile(pdb.topology, nvt_state.getPositions(), f_pdb, keepIds=True)
        with open(nvt_state_path, 'w', encoding='utf-8') as f_state:
            f_state.write(mm.XmlSerializer.serialize(nvt_state))
        log(f"NVT checkpoint saved to {nvt_pdb_path} and {nvt_state_path}.")
        
        # Clean up NVT simulation context to free GPU resources for NPT
        nvt_final_positions = nvt_state.getPositions()
        nvt_final_velocities = nvt_state.getVelocities()
        nvt_final_box = nvt_state.getPeriodicBoxVectors()
        del sim_nvt
        del integrator_nvt
        del system_nvt
        
        # ----------------------------------------------------------------------
        # PART II: NPT DENSITY EQUILIBRATION & RESTRAINT RELEASE
        # ----------------------------------------------------------------------
        log("-" * 80)
        log("PART II: NPT Density Equilibration & Restraint Release")
        log("-" * 80)
        
        # Load fresh system for NPT
        with open(xml_path, 'r', encoding='utf-8') as f:
            system_npt = mm.XmlSerializer.deserialize(f.read())
            
        # Re-add CustomExternalForce with parameter k
        restraint_npt = mm.CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
        initial_npt_k = npt_stages[0]['restraint_k_kcal_per_mol_a2'] * KCAL_PER_MOL_A2_TO_KJ_PER_MOL_NM2
        restraint_npt.addGlobalParameter("k", initial_npt_k)
        restraint_npt.addPerParticleParameter("x0")
        restraint_npt.addPerParticleParameter("y0")
        restraint_npt.addPerParticleParameter("z0")
        
        # Reference coordinates remain the Phase 4B minimized positions
        for idx in rna_heavy_indices:
            pos_nm = pdb.positions[idx].value_in_unit(unit.nanometer)
            restraint_npt.addParticle(idx, [pos_nm[0], pos_nm[1], pos_nm[2]])
        system_npt.addForce(restraint_npt)
        
        # Add Monte Carlo Barostat
        barostat = mm.MonteCarloBarostat(
            npt_press_bar * unit.bar,
            npt_temp_k * unit.kelvin,
            barostat_freq
        )
        barostat.setRandomNumberSeed(random_seed + 1)
        system_npt.addForce(barostat)
        
        # NPT Integrator
        integrator_npt = mm.LangevinMiddleIntegrator(
            npt_temp_k * unit.kelvin,
            friction_ps / unit.picosecond,
            timestep_fs * unit.femtosecond
        )
        integrator_npt.setRandomNumberSeed(random_seed + 2)
        
        sim_npt = app.Simulation(pdb.topology, system_npt, integrator_npt, platform)
        sim_npt.context.setPositions(nvt_final_positions)
        sim_npt.context.setVelocities(nvt_final_velocities)
        sim_npt.context.setPeriodicBoxVectors(*nvt_final_box)
        
        # Prepare NPT CSV
        npt_csv_file = open(npt_csv_path, 'w', encoding='utf-8')
        npt_csv_file.write("stage,step,time_ps,temp_k,pressure_bar,pe_kj_mol,ke_kj_mol,total_energy_kj_mol,volume_nm3,density_g_cm3,restraint_k_kcal_a2\n")
        
        t_npt_start = time.time()
        global_npt_step = 0
        total_npt_sim_ps = nvt_dur_ps
        
        stage_checkpoints = []
        for stage_idx, stage_info in enumerate(npt_stages, 1):
            stage_name = stage_info['name']
            s_dur_ps = stage_info['duration_ps'] * scale_factor
            s_k_kcal = stage_info['restraint_k_kcal_per_mol_a2']
            s_k_kj_nm2 = s_k_kcal * KCAL_PER_MOL_A2_TO_KJ_PER_MOL_NM2
            s_steps = int(round((s_dur_ps * 1000.0) / timestep_fs))
            
            log(f"Beginning NPT {stage_name.upper()} ({stage_idx}/{len(npt_stages)}):")
            log(f"  Duration: {s_dur_ps:.1f} ps ({s_steps} steps) | Restraint k: {s_k_kcal:.2f} kcal/(mol*A^2) ({s_k_kj_nm2:.2f} kJ/(mol*nm^2))")
            
            sim_npt.context.setParameter("k", s_k_kj_nm2)
            
            stage_steps_run = 0
            while stage_steps_run < s_steps:
                steps_chunk = min(report_steps, s_steps - stage_steps_run)
                sim_npt.step(steps_chunk)
                stage_steps_run += steps_chunk
                global_npt_step += steps_chunk
                
                state = sim_npt.context.getState(getEnergy=True, getPositions=True)
                ke = state.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
                pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                inst_T = compute_instantaneous_temperature(ke, system_npt)
                box_vecs = state.getPeriodicBoxVectors().value_in_unit(unit.nanometer)
                vol_nm3, density_g_cm3 = compute_box_density(system_npt, box_vecs)
                time_total_ps = nvt_dur_ps + (global_npt_step * timestep_fs) / 1000.0
                
                npt_csv_file.write(
                    f"{stage_name},{global_npt_step},{time_total_ps:.3f},{inst_T:.2f},{npt_press_bar:.2f},"
                    f"{pe:.3f},{ke:.3f},{pe+ke:.3f},{vol_nm3:.4f},{density_g_cm3:.4f},{s_k_kcal:.2f}\n"
                )
                
                if stage_steps_run % (report_steps * 25) == 0 or stage_steps_run == s_steps:
                    log(f"  [{stage_name.upper()}] Step {stage_steps_run:6d}/{s_steps:6d} (t={time_total_ps:5.1f} ps) | T: {inst_T:5.1f} K | Vol: {vol_nm3:7.3f} nm^3 | Density: {density_g_cm3:6.4f} g/cm^3 | PE: {pe:11.1f} kJ/mol")
                    
            # Save stage checkpoint
            stage_state = sim_npt.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            stage_pdb = os.path.join(npt_dir, f"npt_{stage_name}.pdb")
            stage_xml = os.path.join(npt_dir, f"npt_{stage_name}_state.xml")
            with open(stage_pdb, 'w', encoding='utf-8') as f_pdb:
                app.PDBFile.writeFile(pdb.topology, stage_state.getPositions(), f_pdb, keepIds=True)
            with open(stage_xml, 'w', encoding='utf-8') as f_xml:
                f_xml.write(mm.XmlSerializer.serialize(stage_state))
            stage_checkpoints.append({
                'stage': stage_name,
                'restraint_k_kcal_a2': float(s_k_kcal),
                'duration_ps': float(s_dur_ps),
                'checkpoint_pdb': stage_pdb,
                'checkpoint_xml': stage_xml
            })
            log(f"Checkpoint for {stage_name} written to {stage_pdb}.")
            
        npt_csv_file.close()
        t_npt_elapsed = time.time() - t_npt_start
        log(f"NPT density equilibration completed in {t_npt_elapsed:.2f}s.")
        
        # ----------------------------------------------------------------------
        # PART III: SAVE FINAL EQUILIBRATED ARTIFACTS & VALIDATION
        # ----------------------------------------------------------------------
        final_state = sim_npt.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
        with open(final_pdb_path, 'w', encoding='utf-8') as f_pdb:
            app.PDBFile.writeFile(pdb.topology, final_state.getPositions(), f_pdb, keepIds=True)
        final_pdb_sha = compute_sha256(final_pdb_path)
        with open(final_state_path, 'w', encoding='utf-8') as f_state:
            f_state.write(mm.XmlSerializer.serialize(final_state))
        log(f"Final equilibrated coordinates written to: {final_pdb_path} (SHA-256: {final_pdb_sha})")
        log(f"Final state serialized to: {final_state_path}")
        
        # Automated validation
        log("Executing post-equilibration automated scientific validation...")
        val_results = validate_equilibrated_system(
            final_pdb_path, min_pdb_path, info, npt_csv_path
        )
        for kv, vv in val_results.items():
            if kv not in ['passed', 'errors']:
                log(f"  * {kv}: {vv}")
        if val_results['passed']:
            log("[PASS] All post-equilibration physical, structural, and thermodynamic checks verified!")
        else:
            log(f"[FAIL] Post-equilibration validation failed: {val_results['errors']}")
            raise AssertionError(f"Validation failed for {info['construct_name']}: {val_results['errors']}")
            
        end_time_utc = datetime.datetime.now(datetime.timezone.utc)
        total_walltime = (end_time_utc - start_time_utc).total_seconds()
        
        # Compile metadata
        metadata = {
            'construct': construct_key,
            'construct_id': info['construct_id'],
            'construct_name': info['construct_name'],
            'selected_trna_identifier': info['target_trna'],
            'asl_sequence_5_to_3': info['sequence_5_to_3'],
            'input_files': {
                'system_xml': xml_path,
                'system_xml_sha256': compute_sha256(xml_path),
                'minimized_pdb': min_pdb_path,
                'minimized_pdb_sha256': compute_sha256(min_pdb_path)
            },
            'output_files': {
                'final_equilibrated_pdb': final_pdb_path,
                'final_equilibrated_pdb_sha256': final_pdb_sha,
                'final_state_xml': final_state_path,
                'nvt_pdb': nvt_pdb_path,
                'nvt_csv': nvt_csv_path,
                'npt_csv': npt_csv_path,
                'equilibration_log': log_path
            },
            'stage_checkpoints': stage_checkpoints,
            'openmm_environment': {
                'openmm_version': mm.__version__,
                'platform_used': plat_name,
                'platform_hardware': hw_info,
                'precision_mode': plat_props.get('Precision', 'default'),
                'python_version': sys.version.split()[0]
            },
            'equilibration_protocol': {
                'scale_factor': float(scale_factor),
                'timestep_fs': float(timestep_fs),
                'friction_coeff_ps': float(friction_ps),
                'random_seed': int(random_seed),
                'nvt': {
                    'duration_ps': float(nvt_dur_ps),
                    'temp_start_k': float(nvt_cfg.get('temperature_start_kelvin', 10.0)),
                    'temp_target_k': float(nvt_cfg.get('temperature_target_kelvin', 310.0)),
                    'ramp_duration_ps': float(nvt_ramp_ps),
                    'plateau_duration_ps': float(nvt_plateau_ps),
                    'restraint_k_kcal_a2': float(nvt_k_kcal),
                    'restrained_heavy_atoms': len(rna_heavy_indices)
                },
                'npt': {
                    'duration_ps': float(npt_dur_ps),
                    'target_temp_k': float(npt_temp_k),
                    'target_pressure_bar': float(npt_press_bar),
                    'barostat_frequency_steps': int(barostat_freq),
                    'stages': npt_stages
                },
                'final_restraint_state': 'UNRESTRAINED (k = 0.0 kcal/(mol*A^2))'
            },
            'validation_results': val_results,
            'timing': {
                'start_time_utc': start_time_utc.isoformat(),
                'end_time_utc': end_time_utc.isoformat(),
                'walltime_seconds': float(total_walltime)
            },
            'status': 'PASS',
            'disclaimer': (
                "Equilibration relaxes solvent density and heats the system to target thermodynamic conditions, "
                "but does NOT establish conformational stability, native folding, or biological function."
            )
        }
        
        with open(metadata_yaml_path, 'w', encoding='utf-8') as f_yaml:
            yaml.dump(metadata, f_yaml, sort_keys=False, default_flow_style=False)
        with open(metadata_json_path, 'w', encoding='utf-8') as f_json:
            json.dump(metadata, f_json, indent=2)
            
        log(f"Equilibration metadata exported to {metadata_yaml_path}")
        log(f"Phase 4C Equilibration finished successfully for {info['construct_name']} (Verdict: PASS).\n")
        
    return metadata

def main():
    parser = argparse.ArgumentParser(
        description="Phase 4C: Restrained NVT/NPT Equilibration of WT and CUA Solvated ASL Systems"
    )
    parser.add_argument(
        '--construct', choices=['wt', 'cua', 'all'], default='all',
        help="Construct to equilibrate ('wt', 'cua', or 'all' [default])"
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
        '--scale-factor', type=float, default=1.0,
        help="Scale factor for equilibration durations (default: 1.0 = full 1.2 ns protocol)"
    )
    parser.add_argument(
        '--test-mode', action='store_true',
        help="Run short controlled verification test (scale_factor=0.01, ~12 ps total)"
    )
    parser.add_argument(
        '--validate-only', action='store_true',
        help="Run validation only on existing equilibrated artifacts"
    )
    
    args = parser.parse_args()
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_path = os.path.join(base_dir, args.config)
    
    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
        
    scale_factor = 0.01 if args.test_mode else args.scale_factor
    constructs = ['wt', 'cua'] if args.construct == 'all' else [args.construct]
    
    if args.validate_only:
        print("=" * 80)
        print("PHASE 4C: RUNNING POST-EQUILIBRATION VALIDATION ONLY")
        print("=" * 80)
        all_ok = True
        for c in constructs:
            info = CONSTRUCT_INFO[c]
            p_final = os.path.join(base_dir, 'simulations', info['dir_name'], 'equilibration', 'final_equilibrated.pdb')
            p_min = os.path.join(base_dir, 'simulations', info['dir_name'], 'minimized.pdb')
            p_csv = os.path.join(base_dir, 'simulations', info['dir_name'], 'equilibration', 'npt', 'npt_thermo.csv')
            
            if not os.path.exists(p_final):
                print(f"[ERROR] Missing final_equilibrated.pdb for {c}. Run equilibration first.")
                all_ok = False
                continue
                
            v = validate_equilibrated_system(p_final, p_min, info, p_csv)
            print(f"\nValidation for {info['construct_name']}: {'PASS' if v['passed'] else 'FAIL'}")
            for k, val in v.items():
                if k not in ['passed', 'errors']:
                    print(f"  * {k}: {val}")
            if not v['passed']:
                print(f"  ERRORS: {v['errors']}")
                all_ok = False
        sys.exit(0 if all_ok else 1)
        
    print("=" * 80)
    print("PHASE 4C: SOLVATED ASL NVT/NPT EQUILIBRATION")
    print(f"Constructs: {constructs} | Platform: {args.platform} | Scale Factor: {scale_factor:.4f}")
    print("=" * 80)
    
    results = {}
    for c in constructs:
        results[c] = run_equilibration(
            construct_key=c,
            config=config,
            platform_pref=args.platform,
            scale_factor=scale_factor,
            base_dir=base_dir
        )
        
    if len(constructs) == 2:
        wt_m = results['wt']
        cua_m = results['cua']
        print("\n" + "=" * 80)
        print("PHASE 4C COMPARATIVE EQUILIBRATION SUMMARY (WT vs. CUA)")
        print("=" * 80)
        print(f"{'Metric':<34} | {'Wild-Type (WT)':<20} | {'Suppressor (CUA)':<20}")
        print("-" * 80)
        print(f"{'Sequence (5\'->3\')':<34} | {wt_m['asl_sequence_5_to_3']:<20} | {cua_m['asl_sequence_5_to_3']:<20}")
        print(f"{'Total System Particles':<34} | {wt_m['validation_results']['total_atoms']:<20} | {cua_m['validation_results']['total_atoms']:<20}")
        print(f"{'Platform / Device':<34} | {wt_m['openmm_environment']['platform_used']}: {wt_m['openmm_environment']['platform_hardware'].get('DeviceName', 'CPU'):<14} | {cua_m['openmm_environment']['platform_used']}: {cua_m['openmm_environment']['platform_hardware'].get('DeviceName', 'CPU'):<14}")
        print(f"{'Equilibration Seed':<34} | {wt_m['equilibration_protocol']['random_seed']:<20} | {cua_m['equilibration_protocol']['random_seed']:<20}")
        print(f"{'NVT Duration (ps)':<34} | {wt_m['equilibration_protocol']['nvt']['duration_ps']:<20.1f} | {cua_m['equilibration_protocol']['nvt']['duration_ps']:<20.1f}")
        print(f"{'NPT Duration (ps)':<34} | {wt_m['equilibration_protocol']['npt']['duration_ps']:<20.1f} | {cua_m['equilibration_protocol']['npt']['duration_ps']:<20.1f}")
        print(f"{'Late NPT Mean Temp (K)':<34} | {wt_m['validation_results'].get('npt_late_mean_temperature_k', 0.0):<20.2f} | {cua_m['validation_results'].get('npt_late_mean_temperature_k', 0.0):<20.2f}")
        print(f"{'Late NPT Density (g/cm^3)':<34} | {wt_m['validation_results'].get('npt_late_mean_density_g_cm3', 0.0):<20.4f} | {cua_m['validation_results'].get('npt_late_mean_density_g_cm3', 0.0):<20.4f}")
        print(f"{'Late NPT Mean Volume (nm^3)':<34} | {wt_m['validation_results'].get('npt_late_mean_volume_nm3', 0.0):<20.3f} | {cua_m['validation_results'].get('npt_late_mean_volume_nm3', 0.0):<20.3f}")
        print(f"{'RNA Heavy-Atom RMSD (A)':<34} | {wt_m['validation_results']['rna_heavy_atom_rmsd_angstrom']:<20.4f} | {cua_m['validation_results']['rna_heavy_atom_rmsd_angstrom']:<20.4f}")
        print(f"{'Final Restraint State':<34} | {wt_m['equilibration_protocol']['final_restraint_state']:<20} | {cua_m['equilibration_protocol']['final_restraint_state']:<20}")
        print(f"{'Equilibration Status':<34} | {wt_m['status']:<20} | {cua_m['status']:<20}")
        print("=" * 80)
        print("\n[PHASE 4C COMPLETE] Staged equilibration successfully finished and verified.")
        print("IMPORTANT: No production MD has been executed. Systems are frozen and ready for Phase 5.\n")

if __name__ == '__main__':
    main()

