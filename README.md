# tRNA ASL Structural Stability Study

A reproducible scientific research engineering pipeline investigating the structural integrity and conformational dynamics of an isolated human cytoplasmic tRNA anticodon stem-loop (ASL) upon amber suppressor mutation.

---

## 1. Research Question

> **"When the anticodon of a human cytoplasmic tRNA is mutated to a suppressor anticodon, does the anticodon stem-loop (ASL) retain canonical loop architecture or become destabilized?"**

This study conducts a controlled, replicate-aware molecular dynamics comparison between:
1. **Wild-Type Construct:** An isolated 17-nt human cytoplasmic tRNA ASL with its native anticodon.
   - **Crucial Requirement:** The native anticodon must convert to `CUA` by **exactly two nucleotide substitutions** ($\text{Hamming}(\text{Anticodon}_{\text{WT}}, \text{CUA}) = 2$).
   - Target tRNA identity: **Verified human cytoplasmic tRNA-Asp-GTC-2-1** (GtRNAdb ID: `tRNA-Asp-GTC-2-1`, tRNAscan-SE ID: `chr1.trna92`), anticodon `GUC` $\rightarrow$ `CUA` (Hamming distance = 2). Complete provenance in `data/metadata/selected_trna.yaml` and `docs/11_trna_selection.md`.
2. **Suppressor Construct:** The identical 17-nt ASL containing the `CUA` amber suppressor anticodon at positions 34–36.

---

## 2. Benchmark Simulation Specifications

- **Construct Definition:** 17-nucleotide isolated ASL (5-bp stem, 7-nt loop, unmodified RNA). `[PROMPT MANDATE]`
- **Force Field:** AMBER OL3 (`ff99bsc0chiOL3`). `[PROMPT MANDATE]`
- **Solvent & Ions:** Explicit water, $\approx 0.15\ \text{M}\ \text{KCl}$ (net neutralizing $\text{K}^+$ plus physiological excess). `[PROMPT MANDATE]`
  - **Water Model:** Standard TIP3P is established as the baseline solvent model natively compatible with AMBER OL3. Alternative models (e.g., OPC, TIP4P-EW) remain optional only if accompanied by explicit scientific justification and matched ion parameters.
  - **Parameter Coupling:** Water model and ion parameters form a strictly coupled, internally consistent parameter set (`ionsjc_tip3p` for TIP3P).
- **Thermodynamic Ensemble:** $NPT$ at $T = 310.0\ \text{K}$ ($37^\circ\text{C}$), $P = 1.0\ \text{bar}$. `[PROMPT MANDATE]`
- **Sampling Scheme:** 3 independent replicates per construct, minimum $20.0\ \text{ns}$ production per replicate ($N=6$ simulations, $\ge 120.0\ \text{ns}$ aggregate production). `[PROMPT MANDATE]`

---

## 3. Analysis Framework

1. **Backbone RMSD:** Evaluated separately for the 5-bp stem and 7-nt loop. `[PROMPT MANDATE]`
2. **Per-Residue RMSF:** Backbone root-mean-square fluctuations across residues 1 to 17. `[PROMPT MANDATE]`
3. **Stem Base-Pair Hydrogen-Bond Occupancy:** Time-resolved geometric analysis of the 5 canonical base pairs. `[PROMPT MANDATE]`
4. **Canonical U-Turn Diagnostic:** Distance between U33 N3 (and O2') and phosphate 36 oxygens; $\alpha/\zeta$ backbone reversal angles. `[PROMPT MANDATE]`
5. **3'-Anticodon Stacking:** Centroid distances and inter-planar angles across residues 34–35, 35–36, and 36–37. `[PROMPT MANDATE]`
6. **Additional Justified Metric:** Pending formal justification and selection in Phase 2. Candidates under evaluation include Loop Radius of Gyration ($R_{g,\text{loop}}$), Anticodon Solvent-Accessible Surface Area (SASA), and Stem-Loop Interhelical Tilt Angle. `[PENDING VERIFICATION]`
7. **Replicate-Aware Statistical Synthesis:**
   - **Sample Size Limitation:** $N=3$ replicates per construct provide limited statistical power for traditional significance tests (e.g., a two-sided Mann-Whitney $U$ test cannot achieve $p < 0.05$; minimum achievable $p = 0.10$).
   - **Reporting Methodology:** Primary inferences will be drawn from replicate-level effect sizes (Hedges' $g$ / Cohen's $d$, difference in replicate means/medians) with 95% bootstrap confidence intervals. Trajectory frames are never treated as independent biological replicates (avoiding pseudoreplication).

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

- **Phase 1 (Complete):** Requirements, repository architecture, protocols, and audit corrections finalized.
- **Phase 2 (Complete - Candidate Selected & Frozen):**
  - Public GtRNAdb screening across all 432 human cytoplasmic tRNAs.
  - Selected and verified **human cytoplasmic tRNA-Asp-GTC-2-1** ($\text{Hamming}(\text{GUC}, \text{CUA}) = 2$; 5/5 GC stem).
  - Sequence and metadata frozen in `data/metadata/selected_trna.yaml` and documented in `docs/11_trna_selection.md`.
- **Phase 3 (Pending Next Phase):** System solvation, ion placement, and AMBER OL3 parameterization.
