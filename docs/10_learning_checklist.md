# 10 - Research Engineering Learning Checklist

This checklist synthesizes the core biophysical, molecular dynamics, and statistical concepts required to build, execute, and interpret this study with scientific rigor.

---

## 1. RNA Structural Biology & ASL Architecture

- [ ] **Canonical tRNA Secondary Structure & Sprinzl Numbering**
  - Understand the standard cloverleaf architecture and the universal Sprinzl coordinate system.
  - Locate the 17-nt ASL domain spanning residues 27 to 43.
  - Identify the 5-bp Watson-Crick anticodon stem (pairs 27:43, 28:42, 29:41, 30:40, 31:39).
  - Identify the 7-nt anticodon loop (residues 32, 33, 34, 35, 36, 37, 38).
- [ ] **tRNA Mutation Constraint for Amber Suppression**
  - Understand why the wild-type anticodon must differ from `CUA` by **exactly two nucleotide substitutions** ($\text{Hamming} = 2$).
  - Review how this constraint filters eligible human cytoplasmic tRNA candidates in public databases (GtRNAdb).
- [ ] **The Canonical U-Turn Motif (Quigley & Rich, 1976)**
  - Understand how invariant residue U33 initiates the sharp $180^\circ$ reversal of the phosphodiester backbone.
  - Identify key stabilizing hydrogen bonds: U33 N3-H donor to the non-bridging phosphate oxygens of residue 36 (OP1/OP2), and U33 O2'-H to residue 36/35.
  - Understand characteristic backbone torsions: $\zeta_{33}$ and $\alpha_{34}$ adopting signature gauche/trans conformations.
- [ ] **Anticodon Presentation & 3'-Stacking**
  - Understand how the U-turn projects the anticodon nucleobases (positions 34, 35, 36) outward into solvent.
  - Understand the continuous right-handed stacking of positions 34 $\rightarrow$ 35 $\rightarrow$ 36 $\rightarrow$ 37.
  - Understand why purine 37 (adjacent to anticodon base 36) is critical for reinforcing this stack.
- [ ] **Post-Transcriptional Modifications & Biological Impact**
  - Position 34 (Wobble base): Inosine ($I$), 5-methylcytidine ($m^5C$), thiolated uridines ($mcm^5s^2U$)—roles in restricting ribose pucker (C3'-endo) and enabling/restricting wobble decoding.
  - Position 37: Bulky hypermodifications ($t^6A, ms^2t^6A, i^6A, m^1G$)—crucial for preventing cross-loop base pairing and maintaining open loop architecture.
  - Position 38: Pseudouridine ($\Psi$)—stabilizing the stem-loop closing base pair via an extra hydrogen-bond donor (N1-H).
  - Recognize that simulating *unmodified* RNA removes these evolutionary stabilizing factors, predisposing loops to elevated flexibility in silico.
- [ ] **Ribosome Decoding Center Context**
  - Understand how the ribosomal A-site (conserved 16S/18S rRNA bases) inspects the codon-anticodon mini-helix.
  - Appreciate why isolated ASL structural stability does not equate to cellular nonsense suppression efficiency.

---

## 2. Molecular Dynamics Theory & Simulation Methodology

- [ ] **Empirical Molecular Mechanics Force Fields for RNA**
  - Bonded terms: bond stretching, angle bending, and dihedral torsions (Fourier series).
  - Non-bonded terms: Lennard-Jones 12-6 potential and Coulomb electrostatics.
  - The history of AMBER RNA parameter development: why ff99 suffered from unphysical laddering, how bsc0 fixed $\alpha/\zeta$ torsions, and how OL3 (Zgarbová et al., 2011) corrected the $\chi$ glycosidic dihedral.
- [ ] **Solvation & Water Models for Nucleic Acids**
  - Understand why standard TIP3P serves as the native baseline solvent for AMBER OL3 (`ff99bsc0chiOL3`).
  - Understand the requirement that water model and ion parameters must form a strictly matched, coupled set (e.g., `ionsjc_tip3p` for TIP3P).
  - Learn the trade-offs of alternative 4-point models (OPC, TIP4P-EW) and why using them requires explicit justification and matched ion parameters.
- [ ] **Ionic Strength & Electrostatic Neutralization**
  - RNA as a polyanionic electrolyte (charge $-1e$ per phosphodiester linkage; $-16e$ for a 17-nt hairpin).
  - Net neutralization with $\text{K}^+$ counterions vs. physiological buffer ($\sim 0.15\ \text{M}\ \text{KCl}$).
  - Selection of balanced ion parameters (e.g., Joung-Cheatham) to prevent artificial salt precipitation or unphysical RNA-ion binding.
- [ ] **Thermodynamic Ensembles & Boundary Conditions**
  - $NPT$ ensemble: maintaining physiological temperature ($310\ \text{K}$) and pressure ($1.0\ \text{bar}$).
  - Temperature coupling: Langevin dynamics (friction coefficient $\gamma = 1.0\ \text{ps}^{-1}$).
  - Pressure coupling: Monte Carlo isotropic barostat.
  - Periodic boundary conditions (PBC) and Particle Mesh Ewald (PME) for long-range electrostatics.
- [ ] **Equilibration Strategy for RNA**
  - Restrained minimization to eliminate steric clashes without distorting native A-form geometry.
  - Gradual NVT heating with harmonic position restraints on solute heavy atoms.
  - Stepwise restraint release during NPT density relaxation to avoid shock-induced loop unfolding.

---

## 3. Statistical Analysis & Time-Series Inference

- [ ] **Temporal Autocorrelation in Molecular Dynamics**
  - Understand why trajectory frames sampled at short intervals ($\Delta t = 10\ \text{ps}$) are not statistically independent.
  - Calculating autocorrelation functions $C(\tau)$ and integrated autocorrelation time $\tau_{\text{int}}$.
  - Determining effective independent sample size ($N_{\text{eff}} = N / 2\tau_{\text{int}}$) within a single trajectory.
- [ ] **Block Averaging & Replicate Independence**
  - Grouping autocorrelated time series into stationary blocks larger than $2\tau_{\text{int}}$.
  - Strictly avoiding pseudoreplication: understanding why trajectory frames cannot be pooled across replicates as independent data points.
- [ ] **Small Sample Size Power Constraints ($N=3$)**
  - Understand why traditional two-sided non-parametric tests (such as Mann-Whitney $U$) cannot achieve $p < 0.05$ with $N=3$ vs. $N=3$ (minimum achievable $p = 0.10$).
  - Framing statistical comparisons around replicate-level effect sizes (Hedges' $g$, Cohen's $d$) with 95% bootstrap confidence intervals ($10,000$ iterations).
- [ ] **Scientific Candor & Handling Inconclusive Results**
  - Recognizing when data show overlapping distributions or high replicate divergence.
  - Reporting inconclusive findings without bias or forced narrative.

---

## 4. Software Engineering & Tooling Stack

- [ ] **OpenMM Simulation Framework**
  - Setting up `System`, `Topology`, `Positions`, `Integrator`, and `Simulation` objects.
  - Configuring PME, non-bonded cutoffs, and constraint algorithms (SHAKE/SETTLE).
  - Adding `DCDReporter`, `StateDataReporter`, and `CheckpointReporter`.
  - Managing GPU platform execution and deterministic velocity seeds.
- [ ] **Trajectory Analysis Libraries (MDAnalysis & MDTraj)**
  - Selection algebra: selecting backbone atoms, stem base pairs, and loop residues.
  - Superposition and alignment algorithms (aligning stem before computing loop displacement).
  - Geometric calculations: distance arrays, dihedral angles, hydrogen-bond donor-acceptor geometries, and inertia tensors.
- [ ] **Data Hygiene & Reproducibility Standards**
  - Managing `.gitignore` to prevent raw binary trajectories from polluting git.
  - Exporting structured, lightweight summary datasets (`data/processed/`).
  - Generating publication-quality vector plots (`figures/`).
