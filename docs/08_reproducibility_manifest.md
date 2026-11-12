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

# 2. Recreate Conda environment
conda env create -f environment.yml
conda activate trna-asl-md

# 3. Model 17-nt ASL structures (WT and CUA Suppressor)
python scripts/01_prepare_asl_models.py --config config/simulation_config.yaml

# 4. Solvate and parameterize systems
python scripts/02_setup_system.py --config config/simulation_config.yaml

# 5. Execute energy minimization and equilibration
python scripts/03_equilibrate.py --config config/simulation_config.yaml

# 6. Execute 6 production MD replicates (120 ns total)
python scripts/04_run_production.py --config config/simulation_config.yaml

# 7. Compute trajectory metrics
python scripts/05_analyze_trajectories.py --config config/analysis_config.yaml

# 8. Run replicate-aware statistical synthesis and generate publication figures
python scripts/06_statistical_synthesis.py --config config/analysis_config.yaml
```
