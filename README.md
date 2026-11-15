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
│   ├── 10_learning_checklist.md
│   ├── 18_phase6_trajectory_analysis.md
│   └── 19_phase6_2_compliance_audit.md
└── report/                         # Final formal scientific report and manuscript
```

---

## 5. Scientific Governance & Integrity Principles

- **No Fabricated Data or Premature Conclusions:** No simulation trajectories, results, or biological outcomes are assumed or reported without empirical execution.
- **Exclusively Public Biological Data:** All structural and sequence references originate from open, public databases (GtRNAdb, RCSB PDB, MODOMICS). No proprietary Transloka Bio or KritRNA data is used.
- **Transparent Decision-Making:** Every scientific assumption, parameter choice, and analytical decision is documented in `docs/04_scientific_decisions.md`.
- **Lightweight Repository Hygiene:** Raw trajectories (`.dcd`, `.xtc`, `.nc`) are git-ignored; only lightweight processed tabular summaries, figures, and reproducible scripts are committed.

### 5.1 AI Assistant Disclosure
In full compliance with Transloka Bio assessment guidelines:
- **AI Tool Utilized:** Antigravity (Google DeepMind Gemini-powered agentic coding assistant).
- **Scope of AI Assistance:** AI assistance was employed for boilerplate code generation, Python script orchestration (OpenMM simulation drivers, MDTraj analysis pipelines, NumPy/SciPy statistical testing, Matplotlib/Seaborn plotting), and technical Markdown documentation formatting.
- **Scientific Oversight:** All biological hypotheses, candidate tRNA selection criteria, structural modeling decisions, force field selections, statistical testing architectures, and biophysical interpretations were independently designed, audited, and verified by the author. No automated claims were accepted without empirical execution.

---

## 6. Current Project Status

- **Phase 1 (Complete):** Requirements, repository architecture, protocols, and audit corrections finalized.
- **Phase 2 (Complete - Candidate Selected, Frozen & Reproducibility Audited):**
  - Public GtRNAdb screening across all 432 human cytoplasmic tRNAs.
  - Selected and verified **human cytoplasmic tRNA-Asp-GTC-2-1** ($\text{Hamming}(\text{GUC}, \text{CUA}) = 2$; 5/5 GC stem).
  - Multi-tier ranking formalized and frozen in `data/metadata/selected_trna.yaml` and `docs/11_trna_selection.md`.
- **Phase 3 & 3.1 (Complete - Structural Modeling & Hardening):**
  - 3D starting coordinates constructed from PDB 6UGG (Chain A, 1.95 Å) via rigid-body glycosidic base grafting.
  - 14/14 automated structural validation checks verified across Tiers 1–4 (`structures/wt/initial.pdb`, `structures/cua/initial.pdb`).
- **Phase 4A (Complete - Explicit-Solvent System Preparation):**
  - Cubic periodic box ($55.00\ \text{Å}$) with $> 11.2\ \text{Å}$ solute clearance.
  - Neutralized with $16\ \text{K}^+$ counterions; $150\ \text{mM}$ physiological KCl added ($14\ \text{K}^+/\text{Cl}^-$ pairs).
  - OpenMM systems serialized in `simulations/wt/system/` and `simulations/cua/system/`.
- **Phase 4B (Complete - Energy Minimization):**
  - Two-stage energy minimization (Stage 1 restrained solute $k = 5.0\ \text{kcal/mol/Å}^2$; Stage 2 unrestrained full system).
  - Large potential energy reduction ($\Delta E = -168,778\ \text{kJ/mol}$ WT, $-170,644\ \text{kJ/mol}$ CUA) with preserved RNA fold (RMSD $< 0.9\ \text{Å}$) and zero steric clashes.
  - Minimized coordinates frozen in `simulations/wt/minimized.pdb` and `simulations/cua/minimized.pdb`.
- **Phase 4C (Complete - Restrained NVT/NPT Equilibration & Restraint Release):**
  - Staged 1.2 ns equilibration protocol executed for both WT and CUA (200 ps NVT heating 10 to 310 K, followed by 1000 ps NPT density relaxation across 4 restraint release stages: $k = 2.0 \rightarrow 1.0 \rightarrow 0.2 \rightarrow 0.0\ \text{kcal/mol/Å}^2$).
  - Systems relaxed to stationary target density ($1.0187\ \text{g/cm}^3$ WT, $1.0193\ \text{g/cm}^3$ CUA) and physiological temperature ($310.69\ \text{K}$ WT, $310.35\ \text{K}$ CUA).
  - Native stem fold preserved ($< 0.93\ \text{Å}$ stem RMSD) and zero atomic clashes ($> 2.65\ \text{Å}$ non-adjacent contacts).
  - Final unrestrained coordinates and OpenMM states frozen in `simulations/wt/equilibration/` and `simulations/cua/equilibration/`.
  - Comprehensive report in `docs/16_phase4c_equilibration.md`.
- **Phase 5 (Complete - Production Molecular Dynamics & Run-Level QC):**
  - Six independent production MD simulations executed (3 × WT, 3 × CUA) at $310.0\ \text{K}$ and $1.0\ \text{bar}$ in explicit solvent.
  - Replicates simulated to 20.0 ns each ($120.0\ \text{ns}$ aggregate sampling = $60,000,000$ integration steps = $12,000$ trajectory frames).
  - 100% run-level QC passed across all 6 replicates (atom counts conserved, strictly finite energies, stable physiological temperatures $310.3\text{--}310.4\ \text{K}$, densities $1.018\text{--}1.019\ \text{g/cm}^3$, and readable DCD trajectories).
  - Comprehensive documentation in `docs/17_phase5_production.md`.
- **Phase 6 & 6.1 (Complete - Trajectory Analysis, Structural Comparison & Scientific Audit):**
  - Trajectory-based structural analysis executed across all six production replicates ($120.0\ \text{ns}$ aggregate sampling; $N=3$ WT vs $N=3$ CUA).
  - Evaluated Backbone RMSD (Full ASL, 5-bp stem, loop global, loop internal), per-residue RMSF (positions 27–43), stem base-pair hydrogen-bond persistence (intact-state occupancies $\ge 98\%$ across all 5 monitored pairs), canonical loop U-turn diagnostic interaction (U33(N3) $\cdots$ P36 phosphate), 3'-anticodon base stacking geometry, and nucleotide 34 ribose sugar pucker (pseudorotation phase $P$ and North/East/South states).
  - Replicate-level statistics ($N=3$ vs $N=3$, no frame pseudoreplication) calculated with Welch's $t$-test, Cohen's $d$ effect sizes, and Benjamini-Hochberg FDR $q$-values.
  - Phase 6.1 scientific audit completed: verified U-turn chemistry, stem H-bond language, stacking geometry, Altona-Sundaralingam sugar pucker calculations, and replicate consistency across all six primary observables (`analysis/tables/replicate_consistency_audit.csv`).
  - Sampling sufficiency evaluated across four contiguous 5-ns time blocks: the six 20-ns simulations were numerically stable and provided consistent evidence for substantial WT–CUA differences in several local structural observables (consistent with local loop destabilization/reorganization); however, selected loop observables underwent time-dependent transitions within individual replicates (e.g., U-turn breathing in `wt_rep3`, unstacking in `cua_rep2`/`cua_rep3`), indicating that $3 \times 20\ \text{ns}$ does not establish exhaustive equilibrium sampling of the loop conformational landscape.
  - Six publication-quality figures (`figures/figure1_rmsd.png` through `figures/figure6_sugar_pucker.png`) and machine-readable tables in `analysis/tables/`.
  - Comprehensive documentation in `docs/18_phase6_trajectory_analysis.md`.
- **Phase 6.2 (Complete - Take-Home Requirements Compliance Audit):**
  - Item-by-item compliance audit completed against all original Transloka Bio take-home specifications (`docs/19_phase6_2_compliance_audit.md`).
  - Verified natural modification profiles at positions 34 (Q34/manQ34), 37 ($m^1A37$), and 38 ($m^5C38$) and documented biophysical consequences of their omission.
  - Documented equilibration/production cutoff separation (100% of the 1.2 ns equilibration frames excluded from production analysis).
  - Formulated biophysical justification for nucleotide 34 sugar pucker; stated small-$N$ statistical limitations and sampling requirements for equilibrium claims.
  - Articulated the translational boundary (why ASL MD cannot predict cellular suppression) and specified the physical cell-based dual-luciferase reporter assay.
  - Confirmed AI assistant tool disclosure, path portability, and public repository reproducibility.
- **Phase 7 (Complete - Final Scientific Report & Synthesis):**
  - Prepared submission-ready publication report at [report/report.md](file:///d:/trna-asl-md-stability/report/report.md) ($\le 4$ pages main text; compiled PDF at [report/report.pdf](file:///d:/trna-asl-md-stability/report/report.pdf)).
  - Directly answered core scientific question: CUA suppressor mutation induces substantial local loop reorganization (near-total loss of canonical U-turn, unstacking of 34–35, elevation of $O4'\text{-endo}$ ribose pucker) while the 5-bp helical stem remains locked and persistent ($\ge 98.5\%$ base-pairing).
  - Explicitly articulated statistical hierarchy ($N=3$ independent units, Welch's $t$, Cohen's $d$, FDR correction, 100% directional consistency) and sampling boundaries ($3 \times 20\ \text{ns}$ captures kinetic transitions rather than exhaustive equilibrium thermodynamic ensembles).
  - Addressed biochemical impact of omitted native modifications (Q34/manQ34, $m^1A37$, $m^5C38$).
  - Detailed translational boundary separating isolated ASL simulations from whole-tRNA charging and cellular suppression, and specified a dual-luciferase reporter assay in HEK293T cells.
  - Transparent AI tool disclosure and full open-science reproducibility statement.
