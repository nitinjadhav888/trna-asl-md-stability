# tRNA ASL Structural Stability Study

A reproducible scientific research engineering pipeline investigating the structural integrity and conformational dynamics of an isolated human cytoplasmic tRNA anticodon stem-loop (ASL) upon amber suppressor mutation.

---

## 1. Research Question

> **"When the anticodon of a human cytoplasmic tRNA is mutated to a suppressor anticodon, does the anticodon stem-loop (ASL) retain canonical loop architecture or become destabilized?"**

This study conducts a controlled, replicate-aware molecular dynamics comparison between:
1. **Wild-Type Construct:** An isolated 17-nt human cytoplasmic tRNA ASL with its native anticodon.
2. **Suppressor Construct:** The identical 17-nt ASL containing the `CUA` amber suppressor anticodon.

---

## 2. Benchmark Simulation Specifications

- **Construct Definition:** 17-nucleotide isolated ASL (5-bp stem, 7-nt loop, unmodified RNA).
- **Force Field:** AMBER OL3 (`ff99bsc0chiOL3`).
- **Solvent & Ions:** Explicit water, $\approx 0.15\ \text{M}\ \text{KCl}$ (net neutralizing $\text{K}^+$ plus physiological excess).
- **Thermodynamic Ensemble:** $NPT$ at $T = 310.0\ \text{K}$ ($37^\circ\text{C}$), $P = 1.0\ \text{bar}$.
- **Sampling Scheme:** 3 independent replicates per construct, minimum $20.0\ \text{ns}$ production per replicate ($N=6$ simulations, $\ge 120.0\ \text{ns}$ aggregate production).

---

## 3. Analysis Framework

1. **Backbone RMSD:** Evaluated separately for the 5-bp stem and 7-nt loop.
2. **Per-Residue RMSF:** Backbone root-mean-square fluctuations across residues 1 to 17.
3. **Stem Base-Pair Hydrogen-Bond Occupancy:** Time-resolved geometric analysis of the 5 canonical base pairs.
4. **Canonical U-Turn Diagnostic:** Distance between U33 N3 (and O2') and phosphate 36 oxygens; $\alpha/\zeta$ backbone reversal angles.
5. **3'-Anticodon Stacking:** Centroid distances and inter-planar angles across residues 34–35, 35–36, and 36–37.
6. **Additional Justified Metric:** Loop Radius of Gyration ($R_{g,\text{loop}}$) to quantify global loop expansion or compaction.
7. **Replicate-Aware Statistical Synthesis:** Autocorrelation time estimation ($\tau_{\text{int}}$), effective sample size calculation ($N_{\text{eff}}$), block averaging, bootstrap confidence intervals, and non-parametric hypothesis testing between replicate distributions.

---

## 4. Repository Structure

```text
trna-asl-md-stability/
├── README.md                       # Project overview and workflow documentation
├── LICENSE                         # Open-source MIT license
├── .gitignore                      # Prevents committing large trajectories and binary outputs
├── environment.yml                 # Conda environment specification
├── requirements.txt                # Python dependencies
├── config/
│   ├── simulation_config.yaml      # Simulation parameters and deterministic seeds
│   └── analysis_config.yaml        # Metrics parameters, selections, and statistical thresholds
├── data/
│   ├── raw/                        # Public source sequence/PDB files
│   ├── processed/                  # Versioned, lightweight time-series metric tables
│   └── metadata/                   # Run provenance, hardware logs, and checksums
├── structures/                     # Prepared 17-nt ASL PDB structures
├── simulations/                    # Simulation inputs, topology files, and state checkpoints
├── scripts/                        # Sequential pipeline driver scripts (Phases 2-6)
├── analysis/                       # Modular Python analysis library
├── figures/                        # Publication-quality figures and plots
├── notebooks/                      # Exploratory and audit Jupyter notebooks
├── docs/                           # Comprehensive research and governance documentation
│   ├── 01_project_specification.md
│   ├── 02_requirements_matrix.md
│   ├── 03_architecture.md
│   ├── 04_scientific_decisions.md
│   ├── 05_simulation_protocol.md
│   ├── 06_analysis_plan.md
│   ├── 07_limitations.md
│   ├── 08_reproducibility_manifest.md
│   ├── 09_risk_register.md
│   └── 10_learning_checklist.md
└── report/                         # Final formal scientific report and manuscript
```

---

## 5. Scientific Governance & Integrity Principles

- **No Fabricated Data or Premature Conclusions:** No simulation trajectories, results, or biological outcomes are assumed or reported without empirical execution.
- **Exclusively Public Biological Data:** All structural and sequence references originate from open, public databases (GtRNAdb, RCSB PDB, MODOMICS). No proprietary Transloka Bio or KritRNA data is used.
- **Transparent Decision-Making:** Every scientific assumption, parameter choice, and analytical decision is documented in `docs/04_scientific_decisions.md`.
- **Lightweight Repository Hygiene:** Raw trajectories (`.dcd`, `.xtc`, `.nc`) are git-ignored; only lightweight processed tabular summaries and reproducible scripts are committed.

---

## 6. Current Project Status

- **Phase 1 (Complete):**
  - Requirements inspection and traceability matrix constructed.
  - Repository structure, `.gitignore`, license, and environment configurations initialized.
  - Scientific protocols, architecture, risk register, learning checklist, and configuration placeholders established.
- **Phase 2 (Pending User Approval):**
  - Candidate human cytoplasmic tRNA selection.
  - Structural modeling of 17-nt ASL (WT and CUA suppressor).
