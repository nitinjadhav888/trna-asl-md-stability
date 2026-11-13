# 09 - Project Risk Register

This document identifies, assesses, and outlines mitigation strategies for scientific, computational, reproducibility, and interpretation risks associated with the tRNA ASL structural stability study.

---

## Risk Assessment Matrix Summary

| Severity \ Likelihood | Low | Medium | High |
| :--- | :--- | :--- | :--- |
| **Critical** | — | **RSK-SCI-02** (Timescale / Kinetic Trapping) | — |
| **High** | **RSK-COMP-01** (Numerical / Platform Divergence) | **RSK-SCI-01** (Force Field / Solvent Coupling)<br>**RSK-SCI-04** (Template Bias & Lattice Contacts)<br>**RSK-INT-02** (Extrapolating to Cellular Suppression) | **RSK-INT-01** (Over-Interpreting Small Sample Shifts)<br>**RSK-SCI-06** (Static vs Dynamic Conflation) |
| **Medium** | **RSK-REPR-02** (Library Version Drift)<br>**RSK-SCI-05** (Rigid-Body Grafting Strains) | **RSK-COMP-02** (Simulation Interruptions)<br>**RSK-REPR-01** (Trajectory Storage Bloat) | **RSK-SCI-03** (Terminal Stem Fraying Artifact) |
| **Low** | — | — | — |

---

## Detailed Risk Descriptions & Mitigation Protocols

### 1. Scientific Risks

#### RSK-SCI-01: Force Field & Solvent Parameter Inconsistency
- **Description:** Empirical RNA force fields (AMBER OL3) require matched solvent and ion parameters. Introducing mismatched water models or ion sets can cause artificial RNA-ion clustering, spurious base pairing, or unphysical loop compaction.
- **Likelihood:** Medium | **Impact:** High
- **Mitigation:**
  - Standardize on AMBER OL3 (`ff99bsc0chiOL3`) paired natively with standard TIP3P explicit water as the benchmark baseline.
  - Enforce strictly matched monovalent ion parameters (e.g., Joung-Cheatham `ionsjc_tip3p` for TIP3P).
  - Require documented scientific justification and matched parameters for any alternative solvent model.

#### RSK-SCI-02: Timescale Insufficiency & Kinetic Trapping
- **Description:** A simulation duration of 20 ns per replicate is insufficient to observe full thermodynamic folding/unfolding equilibria ($\mu\text{s}-\text{ms}$ process). The system may remain trapped in its initial starting basin regardless of the true global free energy minimum.
- **Likelihood:** Medium | **Impact:** Critical
- **Mitigation:**
  - Execute 3 independent replicates per construct initialized with different random velocity seeds to probe phase space divergence.
  - Formally assess whether 3 × 20 ns is sufficient in the project report (REQ-20) and explicitly characterize any lack of transition as kinetically limited sampling.

#### RSK-SCI-03: Artificial End-Fraying of the 17-nt Stem Termini
- **Description:** In an isolated 17-nt hairpin, terminal base pair 27:43 lacks flanking base stacking and may fray artificially due to end-solvent exposure or charged terminal groups.
- **Likelihood:** High | **Impact:** Medium
- **Mitigation:**
  - Employ neutral 5'-OH and 3'-OH termini (SDR-09 baseline) to avoid placing repulsive $-2e$ charges on residue 27.
  - Monitor H-bond occupancy for each individual base pair (27:43 through 31:39) separately rather than lumping stem stability into a single scalar average.

#### RSK-SCI-04: Structural Template Bias & Crystalline Lattice Contact Artifacts
- **Description:** Deriving the starting ASL coordinates from crystallographic template PDB 6UGG carries the risk that loop conformations observed in the crystal are stabilized by crystal lattice contacts rather than solution equilibrium thermodynamics. In 6UGG, Chain B adopted an open, non-canonical loop ($13.36\ \text{Å}$ U-turn distance) due to crystal packing, whereas Chain A retained the canonical closed U-turn ($3.52\ \text{Å}$). Selecting Chain A provides the canonical starting baseline, but MD relaxation at 310 K in solvent is required to determine whether this state remains stable without crystalline packing.
- **Likelihood:** Medium | **Impact:** High
- **Mitigation:**
  - Initialize multiple independent replicates ($N=3$) with distinct random seeds to avoid single-trajectory trapping.
  - Continuously track the $U33\ \text{N3} \cdots 36\ \text{OP2}$ distance time series throughout simulation trajectories.
  - Explicitly document that starting from Chain A evaluates the *retention* or *loss* of the canonical fold, rather than assuming it is universally pre-determined.

#### RSK-SCI-05: Rigid-Body Mutation Strains & Local Steric Relic Artifacts
- **Description:** Mutating nucleobases via rigid-body frame grafting ($G34 \rightarrow C34$ and $C36 \rightarrow A36$) places idealized geometry without physical relaxation. While zero hard steric clashes are present ($>3.27\ \text{Å}$ to nearest neighbor), local van der Waals and electrostatic non-bonded interactions are not in their local potential energy minimum prior to Phase 4 minimization.
- **Likelihood:** Low | **Impact:** Medium
- **Mitigation:**
  - Verified zero hard steric clashes ($>1.80\ \text{Å}$ threshold, actual minimum $>3.27\ \text{Å}$).
  - Mandatory multi-stage steepest-descent / L-BFGS energy minimization with restrained solute in Phase 4 before MD production.

#### RSK-SCI-06: Static Starting Geometry vs. Dynamic Stability Conflation
- **Description:** Interpreting initial structure validation (Tier 1–4 checks) as evidence that the CUA suppressor ASL is "stable" or "functional".
- **Likelihood:** High | **Impact:** High
- **Mitigation:**
  - Explicit separation of concepts across all documentation: static validation establishes initial starting sanity; dynamic stability is evaluated by Phase 5 MD; biological suppression efficacy cannot be determined by isolated ASL MD alone.

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

#### RSK-INT-01: Over-Interpreting Small Sample Size Inferences ($N=3$)
- **Description:** Asserting statistically significant differences between WT and Suppressor constructs based on $N=3$ replicates, or committing pseudoreplication by treating autocorrelated trajectory frames as independent data points.
- **Likelihood:** High | **Impact:** High
- **Mitigation:**
  - Explicitly acknowledge that with $N=3$ per group, traditional two-sided non-parametric tests cannot achieve $p < 0.05$ (minimum achievable $p = 0.10$).
  - Primary reporting will rely on standardized effect sizes (Hedges' $g$, Cohen's $d$) with 95% bootstrap confidence intervals.
  - Never pool trajectory frames across replicates as independent observations.
  - If differences show high replicate variance or overlap, report results transparently as inconclusive (REQ-21).

#### RSK-INT-02: Conflating ASL Structural Stability with Cellular Suppression Efficiency
- **Description:** Concluding that because an ASL is structurally stable in MD, the suppressor tRNA will exhibit high readthrough efficiency in human cells.
- **Likelihood:** Medium | **Impact:** High
- **Mitigation:**
  - Address this limitation explicitly in the project specification and final report (REQ-22).
  - Detail the essential cellular processes (aaRS aminoacylation, eEF1A binding, ribosome decoding center discrimination, eRF1 competition) that cannot be modeled by an isolated ASL in solvent.
