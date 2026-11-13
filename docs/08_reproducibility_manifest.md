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
