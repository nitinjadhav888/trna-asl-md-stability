# 08 - Reproducibility & Provenance Manifest

This manifest documents the reproducibility protocols, versioning standards, deterministic seed controls, data provenance policies, and verification steps required to ensure this research pipeline can be audited and reproduced by independent researchers.

---

## 1. Computational Environment & Dependency Pinning

To guarantee software reproducibility across platforms:
- **Conda Environment Specification:** Version-pinned packages documented in `environment.yml`.
- **Python Packaging:** Secondary dependency list maintained in `requirements.txt`.
- **Environment Provenance Capture:** The pipeline records an automated environment snapshot (`data/metadata/environment_info.json`) prior to simulation execution, capturing:
  - Operating System name and kernel release
  - Python interpreter build and architecture
  - OpenMM / AmberTools / MDAnalysis package versions
  - CUDA / OpenCL driver versions and GPU device hardware identifiers

---

## 2. Pseudo-Random Number Generation & Seed Control

Molecular dynamics simulations require stochastic velocity initialization. To enable exact algorithmic replay:
- Maxwell-Boltzmann initial velocity generation uses deterministic integer seeds explicitly configured in `config/simulation_config.yaml`:
  - Wild-Type Replicate 1: `seed = 1001`
  - Wild-Type Replicate 2: `seed = 1002`
  - Wild-Type Replicate 3: `seed = 1003`
  - Suppressor Replicate 1: `seed = 2001`
  - Suppressor Replicate 2: `seed = 2002`
  - Suppressor Replicate 3: `seed = 2003`
- Langevin friction forces and Monte Carlo barostat volume exchange attempts are seeded deterministically per trajectory.

---

## 3. Data Integrity, Provenance & Governance

### 3.1 Public-Only Biological Data Sources
In strict adherence to governance constraints:
- All tRNA sequences, secondary structures, and 3D templates are retrieved solely from open, public biological repositories:
  - **GtRNAdb:** Genomic tRNA Database (Chan & Lowe, UC Santa Cruz)
  - **RCSB Protein Data Bank (PDB):** Experimentally determined 3D atomic coordinates
  - **MODOMICS:** Database of RNA modification pathways (Boccaletto et al.)
  - **NCBI Nucleotide / GenBank**
- **Strict Compliance:** Zero proprietary data from Transloka Bio, KritRNA, or any undisclosed commercial entity is used.

### 3.2 Checksumming & Provenance Tracking
- Every starting PDB structure, topology, and final processed CSV dataset is hashed using SHA-256 upon generation.
- Checksums and generation timestamps are recorded in `data/metadata/manifest_checksums.json`.

---

## 4. Trajectory Archival & Data Hygiene Policy

- **Exclusion of Binary Trajectories:** In accordance with REQ-26, raw binary trajectory files (`.dcd`, `.xtc`, `.nc`, `.trr`) are excluded from git via `.gitignore` to prevent repository bloat.
- **Lightweight Processed Artifacts:** Time-series arrays (RMSD, RMSF, H-bond occupancies, stacking metrics) are exported as compressed CSV/JSON files into `data/processed/`, enabling full reproduction of all downstream figures and statistical tests without re-running MD simulations.
- **Minimal Checkpoint Retention:** Only initial and final simulation state XMLs (`openmm_state_*.xml`) and PDB checkpoints are archived to document simulation convergence.

---

## 5. Independent Verification Protocol

An independent peer reviewer or scientific evaluator can re-run this pipeline from scratch using the following standard sequence:

```powershell
# 1. Clone repository
git clone https://github.com/example/trna-asl-md-stability.git
cd trna-asl-md-stability

# 2. Recreate Conda / Python environment
conda env create -f environment.yml
conda activate trna-asl-md

# 3. Phase 2: Execute biological data verification & candidate selection
python scripts/01_select_trna.py

# 4. Phase 3: Construct 3D ASL structures and validate stereochemistry
python scripts/02_build_asl_structures.py
python scripts/02_validate_asl_structure.py
python scripts/03_visualize_asl_structures.py

# 5. Phase 4: Solvate and parameterize systems (Pending Phase 4)
# python scripts/04_setup_system.py --config config/simulation_config.yaml

# 6. Phase 4: Execute energy minimization and equilibration (Pending Phase 4)
# python scripts/05_equilibrate.py --config config/simulation_config.yaml

# 7. Phase 5: Execute 6 production MD replicates (120 ns total, Pending Phase 5)
# python scripts/06_run_production.py --config config/simulation_config.yaml

# 8. Phase 5: Compute trajectory metrics (Pending Phase 5)
# python scripts/07_analyze_trajectories.py --config config/analysis_config.yaml

# 9. Phase 6: Run replicate-aware statistical synthesis and generate figures
# python scripts/08_statistical_synthesis.py --config config/analysis_config.yaml
```

---

## 6. Phase 3 3D Structural Modeling Provenance & Checksums

- **Date of Execution:** November 13, 2026
- **Modeling Strategy:** Experimental crystallographic template-based homology modeling with rigid-body glycosidic base grafting (`docs/12_structure_modeling_method.md`).
- **Parent Crystallographic Template:** PDB `6UGG` (Chain A, residues 27–43, 1.95 Å resolution unmodified *E. coli* tRNA-Asp, Chan CW, Chetnani B, Mondragón A, *RNA* 2020, 26:278–291, DOI: [10.1261/rna.073478.119](https://doi.org/10.1261/rna.073478.119), PMID: 31848215). Chain A contains the canonical folded U-turn ($d(U33\ \text{N3} \cdots C36\ \text{OP2}) = 3.52\ \text{Å}$), whereas Chain B exhibits an extended non-canonical loop ($13.36\ \text{Å}$) stabilized by crystal lattice contacts.
- **Software Dependencies:**
  - Python: 3.11+
  - BioPython: 1.86
  - NumPy: 2.2.6
  - Matplotlib: 3.10.0
  - OpenMM: 8.6.1 (Force field template verification with `amber14/RNA.OL3.xml` and `amber14/tip3p.xml`)
- **Deterministic Scripts (Phase 3.1 Hardened):**
  - Builder: [scripts/02_build_asl_structures.py](file:///d:/trna-asl-md-stability/scripts/02_build_asl_structures.py) (SHA-256: `a202761b4b6c01c3b229d94943f0d1fce31051a40132300c6dbcfc829ce05019`)
  - Validator: [scripts/02_validate_asl_structure.py](file:///d:/trna-asl-md-stability/scripts/02_validate_asl_structure.py) (SHA-256: `6015a76c30f45bb55743a4b862c82a6620a48a75366b2f08d396a57791dd3507`)
  - Visualizer: [scripts/03_visualize_asl_structures.py](file:///d:/trna-asl-md-stability/scripts/03_visualize_asl_structures.py) (SHA-256: `d28041e375e84e8986bdf087524a2fba339e4083c2d51bd2ff8cce373e50db43`)
- **Artifact Manifest & SHA-256 Hashes:**
  - Raw Template: `data/raw/6UGG.pdb` (SHA-256: `2148a53b3b9e82ed1df163b1be789ab1fb4fcc3b3f461b38905232f4e98ac919`)
  - WT Structure: `structures/wt/initial.pdb` (SHA-256: `61c8363ca99475145ebb266166007e81e32f5bc07aba77ac04906a335b3f291c`)
  - CUA Structure: `structures/cua/initial.pdb` (SHA-256: `277df2aa69ce4240e2434ccfd82dfd387732f10a4054700f44b2dbfe09e31069`)
  - WT Figure: `figures/wt_asl_initial.png` (SHA-256: `8fe2ae10c27a92f3d993a16b145ab51c1f27b37c45e5e9d9c3cbcb1c975a2dfb`)
  - CUA Figure: `figures/cua_asl_initial.png` (SHA-256: `8c9ac1755e5ece3723f0733ce3684f4f87aff2e0a6514929dd71bb219fa8d91b`)
  - Comparison Figure: `figures/asl_structural_comparison.png` (SHA-256: `74810c56691d88b163189c11d24a52c336b9466b24c3eb1936b00656a1dbc75c`)

---

## 7. Phase 4A Solvation & System Preparation Provenance & Checksums

- **Date of Execution:** November 13, 2026
- **System Preparation Strategy:** Explicit TIP3P solvation in an orthogonal cubic periodic boundary ($L = 5.500\ \text{nm} = 55.00\ \text{Å}$) with $16\ \text{K}^+$ neutralizing counterions and $\approx 0.15\ \text{M}$ added physiological KCl ($14\ \text{K}^+ / 14\ \text{Cl}^-$ pairs).
- **Force Fields:**
  - RNA: `amber14/RNA.OL3.xml` (AMBER OL3 / `ff99bsc0chiOL3`)
  - Solvent: `amber14/tip3p.xml` (TIP3P)
  - Ions: `amber14/tip3p.xml` (Joung & Cheatham 2008 parameters for TIP3P, `ionsjc_tip3p`)
- **Software Dependencies:**
  - OpenMM: 8.6.1
  - NumPy: 2.2.6
  - PyYAML: 6.0.2
- **Deterministic Scripts:**
  - System Preparation: [scripts/04_prepare_system.py](file:///d:/trna-asl-md-stability/scripts/04_prepare_system.py) (SHA-256: `bf5668c4bf5e91699e6b18e33d4779cbc42c150038c28cf9d843f84da04cdede`)
  - Scientific Audit: [scripts/audit_phase4a.py](file:///d:/trna-asl-md-stability/scripts/audit_phase4a.py) (SHA-256: `9dcb5ace9de68e932294e53bbdec8119ee325e5e264066cc9ffa0b50e0e8a231`)
- **Artifact Manifest & SHA-256 Hashes:**
  - WT Solvated PDB: `simulations/wt/system/solvated.pdb` (SHA-256: `00b58110d71fcc25be065d744a7402ebff613df0f263663ef4490b0fef37ae1f`)
  - WT System XML: `simulations/wt/system/system.xml` (SHA-256: `d38963ab5ea2adac230bcc552ab1b1dda2260eecf67f41fcf61ceb5a0a8f0b09`)
  - CUA Solvated PDB: `simulations/cua/system/solvated.pdb` (SHA-256: `a3e9c94a497a71245898a0d5c06473c237dbbf63222b7006b31c82d3167813cf`)
  - CUA System XML: `simulations/cua/system/system.xml` (SHA-256: `5e1012f71bde05631b6612f6837c062f7133caf205f3522f3f303fdce7685220`)

---

## 8. Phase 4B Energy Minimization Provenance & Checksums

- **Date of Execution:** November 13, 2026
- **Minimization Strategy:** Two-stage L-BFGS minimization (Stage 1: solute heavy atoms restrained with $k = 5.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$; Stage 2: full system unrestrained, tolerance $10.0\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$).
- **Deterministic Script:** [scripts/05_minimize.py](file:///d:/trna-asl-md-stability/scripts/05_minimize.py) (SHA-256: `a4950196820bf93b96083336b2127b9abdd36ca5ab513f95a99a3c4f86880107`)
- **Artifact Manifest & SHA-256 Hashes:**
  - WT Minimized PDB: `simulations/wt/minimized.pdb` (SHA-256: `471856b2293476ac5bfa24c5b91d43e9ad0476c7aaed8af0bbeacc0ca87efc80`)
  - CUA Minimized PDB: `simulations/cua/minimized.pdb` (SHA-256: `695a1c3beb02b70b84f953b8892f365874dee2756c3df988f5dfa48b752eaf94`)

---

## 9. Phase 4C Restrained NVT/NPT Equilibration Provenance & Checksums

- **Date of Execution:** November 14, 2026
- **Equilibration Strategy:** Staged NVT restrained heating (10 to 310 K, 200 ps, $k = 5.0\ \text{kcal/mol/Å}^2$) followed by 4-stage NPT density equilibration (310 K, 1.0 bar, 1000 ps total, $k = 2.0 \rightarrow 1.0 \rightarrow 0.2 \rightarrow 0.0\ \text{kcal/mol/Å}^2$, 250 ps per stage).
- **OpenMM Platform:** OpenCL (`mixed` precision) on `NVIDIA GeForce MX350`.
- **Deterministic Seeds:**
  - Wild-Type (`wt`): `seed = 4101`
  - Suppressor (`cua`): `seed = 4201`
- **Deterministic Driver Script:** [scripts/06_equilibrate.py](file:///d:/trna-asl-md-stability/scripts/06_equilibrate.py) (SHA-256: `b3248984f04279cdb8084a11825b5d7b4c7714e19c606a16d2b67606cdcd6ddf`)
- **Artifact Manifest & SHA-256 Hashes:**
  - WT Final Equilibrated PDB: `simulations/wt/equilibration/final_equilibrated.pdb` (SHA-256: `84a6f2f690ba7bf013572f7d4a9a91463bec984a5a6703e2ff8f30875cb03409`)
  - WT Final State XML: `simulations/wt/equilibration/final_state.xml` (SHA-256: `0fd3e8b81874636d3683d35f39605ca09c88b846e288ee9c0e8472949622000c`)
  - WT NPT Thermo CSV: `simulations/wt/equilibration/npt/npt_thermo.csv` (SHA-256: `b863aa98651a0b3c66f50c0587d4677764fc8b6dd74c8b201a4e21a24329241b`)
  - CUA Final Equilibrated PDB: `simulations/cua/equilibration/final_equilibrated.pdb` (SHA-256: `e17dd30100bdde7974cac01ef46224283ef65eaab1b683b2adad59766043c65b`)
  - CUA Final State XML: `simulations/cua/equilibration/final_state.xml` (SHA-256: `41b6cfe8ce3cb916aae7242342651df20f9c6ebae768dc9c13ec6d32f7ccb904`)
  - CUA NPT Thermo CSV: `simulations/cua/equilibration/npt/npt_thermo.csv` (SHA-256: `65da5a0f296336a5b6dca74251cb12c9b4e05b5f25a075e7a9b0c034ff66ea9d`)

---

## 10. Phase 5 Independent Production MD Execution & Provenance

- **Date of Execution:** November 14, 2026
- **Protocol:** Isothermal-isobaric ($NPT$) production MD at $T = 310.0\ \text{K}$, $P = 1.0\ \text{bar}$ (Monte Carlo barostat every 25 steps), Langevin Middle integrator ($\gamma = 1.0\ \text{ps}^{-1}$, $\Delta t = 2.0\ \text{fs}$), PME electrostatics ($1.0\ \text{nm}$ cutoff), 120.0 ns aggregate sampling ($6 \times 20.0\ \text{ns}$ = $60,000,000$ steps = $12,000$ DCD frames).
- **Compute Hardware:** NVIDIA Tesla T4 GPU via OpenCL/CUDA (`Precision: mixed`).
- **Deterministic Replicate Seeds:**
  - Wild-Type (`wt`): `wt_rep1` (5101), `wt_rep2` (5102), `wt_rep3` (5103)
  - CUA Suppressor (`cua`): `cua_rep1` (5201), `cua_rep2` (5202), `cua_rep3` (5203)
- **Production Driver Script:** [scripts/07_run_production.py](file:///d:/trna-asl-md-stability/scripts/07_run_production.py)
- **Quality Control Verdict:** 6/6 Replicates PASSED all run-level QC gates (see [docs/17_phase5_production.md](file:///d:/trna-asl-md-stability/docs/17_phase5_production.md)).
- **Repository Hygiene:** Raw trajectory binaries (`*.dcd`) and checkpoints (`checkpoint/`) are excluded via `.gitignore`. Processed thermodynamic logs (`production_thermo.csv`), run logs (`production.log`), final structures (`final_state/`), and metadata are tracked.

---

## 11. Phase 6 Trajectory Analysis & Structural Comparison Provenance

- **Date of Execution:** November 14, 2026
- **Analysis Pipeline & Automation:**
  - Audit: [scripts/08_audit_production_data.py](file:///d:/trna-asl-md-stability/scripts/08_audit_production_data.py)
  - PBC Imaging: [scripts/09_extract_trajectory.py](file:///d:/trna-asl-md-stability/scripts/09_extract_trajectory.py)
  - Backbone RMSD: [scripts/10_rmsd.py](file:///d:/trna-asl-md-stability/scripts/10_rmsd.py)
  - Per-Residue RMSF: [scripts/11_rmsf.py](file:///d:/trna-asl-md-stability/scripts/11_rmsf.py)
  - Stem H-Bond Occupancy: [scripts/12_hbond_analysis.py](file:///d:/trna-asl-md-stability/scripts/12_hbond_analysis.py)
  - Loop U-Turn Diagnostic: [scripts/13_uturn_analysis.py](file:///d:/trna-asl-md-stability/scripts/13_uturn_analysis.py)
  - 3'-Anticodon Base Stacking: [scripts/14_stacking_analysis.py](file:///d:/trna-asl-md-stability/scripts/14_stacking_analysis.py)
  - Sugar Pucker at Nucleotide 34: [scripts/15_additional_metric.py](file:///d:/trna-asl-md-stability/scripts/15_additional_metric.py)
  - Replicate Statistics & Sampling Sufficiency: [scripts/16_statistics.py](file:///d:/trna-asl-md-stability/scripts/16_statistics.py)
  - Publication Figures: [scripts/17_generate_figures.py](file:///d:/trna-asl-md-stability/scripts/17_generate_figures.py)
  - Replicate Consistency & Sampling Convergence Audit: [scripts/18_audit_replicate_consistency.py](file:///d:/trna-asl-md-stability/scripts/18_audit_replicate_consistency.py)
- **Primary Data Artifacts Generated:**
  - Machine-readable summary and time-series tables in [analysis/tables/](file:///d:/trna-asl-md-stability/analysis/tables) (including `replicate_consistency_audit.csv`)
  - Publication figures in [figures/](file:///d:/trna-asl-md-stability/figures) (`figure1_rmsd.png` through `figure6_sugar_pucker.png`)
  - Full scientific analysis document: [docs/18_phase6_trajectory_analysis.md](file:///d:/trna-asl-md-stability/docs/18_phase6_trajectory_analysis.md)
- **Phase Boundary & Scientific Governance Verification:**
  - All comparisons evaluated using $N=3$ WT and $N=3$ CUA independent simulation replicates (Welch's $t$-test, Cohen's $d$, Benjamini-Hochberg FDR).
  - Trajectory frames were not pooled to inflate degrees of freedom.
  - Phase 6 is restricted to structural observables only; no cellular suppression efficiency or biological claims were made.
  - Phase 6.1 scientific audit verified U-turn chemistry, stem H-bond language, stacking geometry, Altona-Sundaralingam sugar pucker calculations, and explicitly documented sampling convergence limitations.

---

## 12. Phase 6.2 Take-Home Requirements Compliance Audit Provenance

- **Date of Execution:** November 15, 2026
- **Audit Scope:** Systematic compliance evaluation against all original Transloka Bio take-home assessment specifications.
- **Audit Documentation:** [docs/19_phase6_2_compliance_audit.md](file:///d:/trna-asl-md-stability/docs/19_phase6_2_compliance_audit.md)
- **Key Verifications:**
  - Native modification profiles audited for positions 34 (Q34/manQ34), 37 ($m^1A37$), and 38 ($m^5C38$).
  - Complete temporal separation between 1.2 ns equilibration and 20.0 ns production trajectories (zero equilibration frames in production analysis).
  - Explicit rationale for position-34 ribose sugar pucker selection documented.
  - Small-$N$ statistical evidence hierarchy and requirements for multi-microsecond equilibrium sampling specified.
  - Translational boundary defined (distinguishing ASL stability, whole-tRNA charging, accommodation, and eRF1 competition) and cell-based dual-luciferase reporter assay in HEK293T cells specified.
  - AI assistant usage fully disclosed; zero proprietary Transloka/KritRNA material present.
- **Phase 7 Readiness:** **PASS** (100% compliant across all 16 prompt requirement areas).

---

## 13. Phase 7 Final Scientific Report and Synthesis Provenance

- **Date of Execution:** November 15, 2026
- **Deliverables:**
  - Final Publication-Ready Report: [report/report.md](file:///d:/trna-asl-md-stability/report/report.md) ($\le 4$ pages main text when rendered; ~2,000 words prose, ~2,600 words with master table).
  - Compiled Report PDF: [report/report.pdf](file:///d:/trna-asl-md-stability/report/report.pdf) (6 pages total including all 6 embedded high-resolution figures).
  - High-Resolution Report Figures: `report/figures/figure1_rmsd.png` through `report/figures/figure6_sugar_pucker.png`.
- **Core Research Question Answered:**
  - Canonical loop architecture is *not* retained in the CUA amber suppressor; the loop undergoes substantial local structural reorganization.
  - The 5-bp helical stem remains completely intact, locked, and persistent ($\ge 98.5\%$ base-pairing, backbone RMSD $< 1.2\ \text{Å}$).
  - Specific localized perturbations: near-total collapse of the U33–P36 U-turn interaction ($87.60 \pm 17.42\%$ in WT vs. $0.03 \pm 0.03\%$ in CUA; distance expands $+2.77\ \text{Å}$), severe uncoupling of 34–35 dinucleotide base stacking ($90.92 \pm 5.93\%$ in WT vs. $19.25 \pm 19.80\%$ in CUA), and a 4-fold elevation in position-34 East ($O4'\text{-endo}$) intermediate ribose sugar pucker ($4.88 \pm 1.70\%$ in WT vs. $19.37 \pm 3.36\%$ in CUA).
- **Statistical Rigor & Sampling Hierarchy:**
  - Inference conducted strictly at replicate level ($N=3$ independent production runs per construct, $120.0\ \text{ns}$ aggregate sampling).
  - Zero frame-level pseudoreplication; Welch's $t$-test ($df \approx 2\text{--}4$), Cohen's $d$ effect sizes, Benjamini-Hochberg FDR correction.
  - 100% directional consistency across all three independent random velocity seeds.
  - Explicit sampling boundary: $3 \times 20\ \text{ns}$ demonstrates consistent local loop remodeling within a 20-ns window, but observed kinetic transitions (`wt_rep3` breathing, `cua_rep2`/`cua_rep3` unstacking) confirm that microsecond-scale enhanced sampling is required for equilibrium thermodynamic claims.
- **Biochemical & Translational Synthesis:**
  - Addressed biochemical consequences of omitted native modifications (Q34/manQ34, $m^1A37$, $m^5C38$).
  - Articulated the translational boundary separating isolated ASL MD from whole-tRNA charging and cellular nonsense suppression.
  - Specified a concrete experimental measurement: cell-based dual-luciferase UAG suppression reporter assay in HEK293T cells.
- **Scientific Integrity & AI Disclosure:**
  - AI assistance (Antigravity / Gemini) disclosed transparently.
  - Zero proprietary Transloka or KritRNA material utilized.
  - All computational pipelines fully documented and reproducible under the MIT License.
- **Final Status:** **PHASE 7 COMPLETE — SUBMISSION READY.**





