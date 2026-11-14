# 17 - Phase 5 Independent Production Molecular Dynamics Execution & QC Report

**Study Title:** tRNA ASL Structural Stability Study  
**Phase:** Phase 5 (Production Molecular Dynamics & Run-Level QC)  
**Status:** Complete — Final Verdict: **PASS (6/6 Replicates Verified)**  
**Aggregate Sampling:** 120.0 ns ($6 \times 20.0\ \text{ns}$ = $60,000,000$ integration steps = $12,000$ trajectory frames)  
**Execution Date:** November 14, 2026  
**Pipeline Driver:** [scripts/07_run_production.py](file:///d:/trna-asl-md-stability/scripts/07_run_production.py)  

---

## 1. Scientific Purpose & Protocol Scope

The objective of Phase 5 is to execute **independent production molecular dynamics simulations** in explicit solvent under isothermal-isobaric ($NPT$) conditions across triplicates for both Wild-Type and CUA Suppressor 17-nt anticodon stem-loops (ASLs).

### Production MD Protocol Specifications
- **Ensemble:** $NPT$ (Isothermal-Isobaric)
- **Target Temperature:** $310.0\ \text{K}$ ($37.0^\circ\text{C}$, physiological human core temperature)
- **Target Pressure:** $1.0\ \text{bar}$ ($1.01325\ \text{bar}$, physiological atmospheric pressure)
- **Thermostat:** OpenMM `LangevinMiddleIntegrator` with collision frequency $\gamma = 1.0\ \text{ps}^{-1}$ and integration timestep $\Delta t = 2.0\ \text{fs}$ (all covalent bonds involving hydrogen constrained with CCMA tolerance $10^{-5}$).
- **Barostat:** OpenMM `MonteCarloBarostat` configured for isotropic periodic box volume scaling attempted every 25 steps ($50\ \text{fs}$).
- **Nonbonded Electrostatics:** Particle Mesh Ewald (PME) with real-space cutoff $1.0\ \text{nm}$ ($10.0\ \text{Å}$) and Ewald error tolerance $0.0005$.
- **Van der Waals:** Lennard-Jones with force-switch starting at $0.9\ \text{nm}$ and cutoff at $1.0\ \text{nm}$.
- **Trajectory Recording:** DCD format recorded every 5,000 steps ($10.0\ \text{ps}$ = 2,000 frames per 20 ns replicate).
- **Thermodynamic Logging:** CSV recording every 5,000 steps ($10.0\ \text{ps}$ = 2,000 records per 20 ns replicate).
- **Periodic State Checkpointing:** Serialized OpenMM State XML and PDB files saved every 500,000 steps ($1.0\ \text{ns}$).
- **Hardware Platform:** NVIDIA Tesla T4 GPU via OpenCL/CUDA (`Precision: mixed`).

### Strict Phase 5 Scope Boundary & Non-Overreach Mandate
- **Mandate:** Phase 5 is strictly limited to **unbiased production trajectory generation, data integrity reconciliation, and run-level quality control (QC)**.
- **Explicit Boundary:** No trajectory post-processing (RMSD time-series, per-residue RMSF profiles, base-pair hydrogen bond occupancies, U-turn geometry distributions, or 3'-base stacking metrics) is reported or inferred within Phase 5. All structural and comparative trajectory analyses are strictly reserved for Phase 6.

---

## 2. Replicate Provenance & Simulation Matrix

All simulations were initiated strictly from the verified, unrestrained Phase 4C final equilibrated coordinate states and Phase 4A system XML definitions:
- **Wild-Type Starting PDB:** [simulations/wt/equilibration/final_equilibrated.pdb](file:///d:/trna-asl-md-stability/simulations/wt/equilibration/final_equilibrated.pdb) (SHA-256: `6fa73db26ba35be196f30ea515d18d8442b5883ef59c5d0b471676916e7aa2e2`)
- **CUA Suppressor Starting PDB:** [simulations/cua/equilibration/final_equilibrated.pdb](file:///d:/trna-asl-md-stability/simulations/cua/equilibration/final_equilibrated.pdb) (SHA-256: `d9ee484511d51a66270b2848c8df37cf2ab479be964893796d11124237fdb9cf`)

### Six-Replicate Production Design

| Construct | Replicate | Seed | Target Duration | Integration Steps | Starting State | Output Directory |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Wild-Type (`wt`)** | Rep 1 | 5101 | 20.0 ns | 10,000,000 | Phase 4C Final Equilibrated | `simulations/wt_rep1` |
| **Wild-Type (`wt`)** | Rep 2 | 5102 | 20.0 ns | 10,000,000 | Phase 4C Final Equilibrated | `simulations/wt_rep2` |
| **Wild-Type (`wt`)** | Rep 3 | 5103 | 20.0 ns | 10,000,000 | Phase 4C Final Equilibrated | `simulations/wt_rep3` |
| **CUA Suppressor (`cua`)** | Rep 1 | 5201 | 20.0 ns | 10,000,000 | Phase 4C Final Equilibrated | `simulations/cua_rep1` |
| **CUA Suppressor (`cua`)** | Rep 2 | 5202 | 20.0 ns | 10,000,000 | Phase 4C Final Equilibrated | `simulations/cua_rep2` |
| **CUA Suppressor (`cua`)** | Rep 3 | 5203 | 20.0 ns | 10,000,000 | Phase 4C Final Equilibrated | `simulations/cua_rep3` |

---

## 3. Run-Level Automated Quality Control Verification

Run-level QC was executed across all six production replicates using `scripts/07_run_production.py --validate-qc`. Every replicate satisfied all eight predetermined QC gates without exception:

| Replicate Tag | Total Atoms | RNA Sequence | Simulated Time | Trajectory Frames | Thermo Records | Mean Temp ($K$) | Mean Density ($g/cm^3$) | NaNs/Infs | Final Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`wt_rep1`** | 15,693 | `CCCGCCUGUCACGCGGG` | 20.0 ns | 2,000 / 2,000 | 2,000 | $310.31 \pm 2.43$ | $1.0187 \pm 0.0041$ | None | **PASS** ✅ |
| **`wt_rep2`** | 15,693 | `CCCGCCUGUCACGCGGG` | 20.0 ns | 2,000 / 2,000 | 2,000 | $310.40 \pm 2.46$ | $1.0184 \pm 0.0039$ | None | **PASS** ✅ |
| **`wt_rep3`** | 15,693 | `CCCGCCUGUCACGCGGG` | 20.0 ns | 2,000 / 2,000 | 2,000 | $310.38 \pm 2.49$ | $1.0185 \pm 0.0040$ | None | **PASS** ✅ |
| **`cua_rep1`** | 15,686 | `CCCGCCUCUAACGCGGG` | 20.0 ns | 2,000 / 2,000 | 2,000 | $310.30 \pm 2.46$ | $1.0184 \pm 0.0039$ | None | **PASS** ✅ |
| **`cua_rep2`** | 15,686 | `CCCGCCUCUAACGCGGG` | 20.0 ns | 2,000 / 2,000 | 2,000 | $310.39 \pm 2.46$ | $1.0184 \pm 0.0039$ | None | **PASS** ✅ |
| **`cua_rep3`** | 15,686 | `CCCGCCUCUAACGCGGG` | 20.0 ns | 2,000 / 2,000 | 2,000 | $310.31 \pm 2.51$ | $1.0183 \pm 0.0039$ | None | **PASS** ✅ |

### Key Run-Level Findings
1. **Atom Count Conservation:** Atom counts remained strictly constant throughout all runs ($15,693$ for WT; $15,686$ for CUA, consistent with the G34C and U36A net mass alteration).
2. **Thermal Stability:** Temperatures across all replicates fluctuated naturally around $310.0\ \text{K}$ with standard deviations $\sigma_T \approx 2.45\text{--}2.51\ \text{K}$, consistent with the canonical Maxwell-Boltzmann distribution for an $N \approx 15,700$ finite particle system in $NPT$.
3. **Physical Box Density:** Bulk densities stabilized at $1.0184\text{--}1.0187\ \text{g/cm}^3$, precisely matching the physical target density for $150\ \text{mM}$ aqueous KCl at $310\ \text{K}$.
4. **Trajectory Continuity:** Trajectories were reconciled and verified for monotonic frame progression ($10.0\ \text{ps}$ intervals, $\Delta t = 2.0\ \text{fs}$) with exactly $2,000$ frames per $20.0\ \text{ns}$ run ($12,000$ frames total).

---

## 4. Sampling Limitations & Explicit Scientific Disclaimers

1. **Finite Timescale (20 ns per replicate):** $20\ \text{ns}$ per replicate ($120\ \text{ns}$ aggregate) is sufficient to observe sub-nanosecond loop fluctuations, local base fraying, and immediate U-turn stability/collapse, but cannot capture rare microsecond-to-millisecond global refolding or complete loop unfolding events.
2. **Classical Empirical Force Field Limitations:** Simulations utilized the AMBER OL3 RNA force field with TIP3P water and Joung-Cheatham monovalent ion parameters. While OL3 is the accepted benchmark for canonical A-form RNA helices, non-canonical loop dynamics can exhibit minor force-field dependent kinetic barriers.
3. **Unmodified Nucleoside Baseline:** All simulations model unmodified standard RNA nucleotides. In native human cells, position 34, 37, and adjacent residues often undergo enzymatic post-transcriptional modifications (e.g., pseudouridylation, queuosine, or methylations) that further modulate loop rigidity.

---

Phase 5 is complete. No Phase 6 or later analysis was performed.
