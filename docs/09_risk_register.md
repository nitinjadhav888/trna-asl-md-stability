# 09 - Project Risk Register

This document identifies, assesses, and outlines mitigation strategies for scientific, computational, reproducibility, and interpretation risks associated with the tRNA ASL structural stability study.

---

## Risk Assessment Matrix Summary

| Severity \ Likelihood | Low | Medium | High |
| :--- | :--- | :--- | :--- |
| **Critical** | — | **RSK-SCI-02** (Timescale/Kinetic Trapping) | — |
| **High** | **RSK-COMP-01** (Numerical/GPU divergence) | **RSK-SCI-01** (Force field loop bias)<br>**RSK-INT-02** (Extrapolating to cellular suppression) | **RSK-INT-01** (Over-interpreting non-significant shifts) |
| **Medium** | **RSK-REPR-02** (Library version drift) | **RSK-COMP-02** (Simulation crash/memory leak)<br>**RSK-REPR-01** (Trajectory storage bloat) | **RSK-SCI-03** (Terminal stem fraying artifact) |
| **Low** | — | — | — |

---

## Detailed Risk Descriptions & Mitigation Protocols

### 1. Scientific Risks

#### RSK-SCI-01: Force Field Bias & Unphysical Loop Compaction
- **Description:** Empirical nucleic acid force fields historically exhibit tendencies toward unphysical base-stacking or non-canonical hydrogen bonding in flexible single-stranded RNA loops, potentially obscuring true physiological dynamics.
- **Likelihood:** Medium | **Impact:** High
- **Mitigation:**
  - Standardize on AMBER OL3 (`ff99bsc0chiOL3`), which was specifically re-parameterized to eliminate non-canonical $\alpha/\zeta$ and $\chi$ rotamer artifacts.
  - Pair OL3 with modern water models (OPC or TIP4P-EW per SDR-06) that reduce artificial hydrophobic compaction relative to TIP3P.
  - Continuously track U-turn hydrogen-bonding and backbone dihedral metrics to flag unphysical conformations early.

#### RSK-SCI-02: Timescale Insufficiency & Kinetic Trapping
- **Description:** A simulation duration of 20 ns per replicate is insufficient to observe full thermodynamic folding/unfolding equilibria ($\mu\text{s}-\text{ms}$ process). The system may remain trapped in its initial starting basin regardless of the true global free energy minimum.
- **Likelihood:** Medium | **Impact:** Critical
- **Mitigation:**
  - Execute 3 independent replicates per construct initialized with different random velocity seeds to probe divergence across phase space.
  - Formally assess whether 3 × 20 ns is sufficient in the project report (REQ-20) and explicitly characterize any lack of transition as kinetically limited sampling.

#### RSK-SCI-03: Artificial End-Fraying of the 17-nt Stem Termini
- **Description:** In an isolated 17-nt hairpin, terminal base pair 27:43 lacks flanking base stacking and may fray artificially due to end-solvent exposure or charged terminal groups.
- **Likelihood:** High | **Impact:** Medium
- **Mitigation:**
  - Employ neutral 5'-OH and 3'-OH termini (SDR-09) to avoid placing repulsive $-2e$ charges on residue 27.
  - Monitor H-bond occupancy for each individual base pair (27:43 through 31:39) separately rather than lumping stem stability into a single scalar average.

---

### 2. Computational & Numerical Risks

#### RSK-COMP-01: Hardware & Numerical Platform Divergence
- **Description:** Running MD on different GPU hardware (NVIDIA CUDA vs. OpenCL vs. CPU) or utilizing mixed-precision floating-point arithmetic can introduce subtle numerical trajectories divergences.
- **Likelihood:** Low | **Impact:** High
- **Mitigation:**
  - Enforce standardized integration precision (e.g., OpenMM `mixed` or `double` precision mode).
  - Record exact hardware devices, driver versions, and compute platforms in metadata logs (`data/metadata/environment_info.json`).

#### RSK-COMP-02: Hardware Interruptions or Simulation Crashes
- **Description:** Hardware faults, power loss, or operating system updates could terminate long production runs prematurely.
- **Likelihood:** Medium | **Impact:** Medium
- **Mitigation:**
  - Write periodic checkpoint state files (`.chk` or `.xml`) every 1.0 ns of simulation time.
  - Structure production scripts to support automatic resumption from the latest valid checkpoint.

---

### 3. Reproducibility Risks

#### RSK-REPR-01: Trajectory Storage Bloat & Git Pollution
- **Description:** Binary MD trajectory files are large (several gigabytes) and committing them to git corrupts repository hygiene and breaches remote repository limits.
- **Likelihood:** High | **Impact:** Medium
- **Mitigation:**
  - Strict `.gitignore` rules actively blocking all binary trajectory and checkpoint extensions (`.dcd`, `.xtc`, `.nc`, `.chk`, `.rst7`).
  - Extract lightweight, processed time-series data tables into `data/processed/` for version-controlled distribution and reporting.

#### RSK-REPR-02: Dependency Version Drift & Environment Incompatibility
- **Description:** Future updates to underlying Python packages (e.g., MDAnalysis, NumPy 2.x breaking changes) could break analysis scripts.
- **Likelihood:** Low | **Impact:** Medium
- **Mitigation:**
  - Provide strict version bounds in `environment.yml` and `requirements.txt` (e.g., `numpy>=1.26.0,<2.0.0`).
  - Containerization / conda environment lockfile support.

---

### 4. Interpretation & Reporting Risks

#### RSK-INT-01: Over-Interpreting Non-Significant Structural Shifts
- **Description:** Observing small numerical differences in RMSD or RMSF between WT and Suppressor and asserting that the suppressor is "destabilized" without statistical significance.
- **Likelihood:** High | **Impact:** High
- **Mitigation:**
  - Mandate replicate-aware statistical hypothesis testing (Mann-Whitney U or permutation testing across replicate means).
  - Explicit rule: If differences do not reach statistical significance ($p < 0.05$ after multiplicity correction) or if replicate variance exceeds between-construct variance, report results as inconclusive (REQ-21).

#### RSK-INT-02: Conflating ASL Structural Stability with Cellular Suppression Efficiency
- **Description:** Concluding that because an ASL is structurally stable in MD, the suppressor tRNA will exhibit high readthrough efficiency in human cells.
- **Likelihood:** Medium | **Impact:** High
- **Mitigation:**
  - Address this limitation explicitly in the project specification and final report (REQ-22).
  - Detail the essential cellular processes (aaRS aminoacylation, eEF1A binding, ribosome decoding center discrimination, eRF1 competition) that cannot be modeled by an isolated ASL in solvent.
