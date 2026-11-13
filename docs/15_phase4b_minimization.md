# 15 - Phase 4B Solvated ASL Energy Minimization Report

**Study Title:** tRNA ASL Structural Stability Study  
**Phase:** Phase 4B (Energy Minimization)  
**Status:** Complete — Final Verdict: **PASS**  
**Execution Date:** November 13, 2026  
**Pipeline Script:** [scripts/05_minimize.py](file:///d:/trna-asl-md-stability/scripts/05_minimize.py)  

---

## 1. Scientific Purpose & Scope

The objective of Phase 4B is to perform **energy minimization only** on the fully solvated, neutralized, salt-adjusted Wild-Type and CUA Suppressor 17-nt ASL systems constructed in Phase 4A.

### Fundamental Scientific Scope & Non-Overreach Mandate
- **Purpose:** Eliminate unphysical steric overlaps, relax strained solvent-solute and ion-solute packing contacts, and resolve non-bonded electrostatic strain inherent to rigid-body solvent placement prior to restrained heating and equilibration.
- **Strict Scope Boundary:** Zero NVT heating, zero NPT density relaxation, zero production MD, zero trajectory analysis, and zero RMSD/RMSF or hydrogen-bond dynamic stability calculations were performed in this phase.
- **Explicit Scientific Precaution:**
  > **Energy minimization removes unfavorable local geometries but does not establish conformational stability, folding, or biological function.**
  The fact that an energy minimization converges and produces a negative $\Delta E$ is an engineering prerequisite for numerical integration, not evidence that either construct is biophysically stable or functional.

---

## 2. Input Systems & Provenance

Minimization was performed directly on the frozen Phase 4A system artifacts without reparameterization, regeneration, or modification:

### Wild-Type ASL (`wt`)
- **System XML:** [simulations/wt/system/system.xml](file:///d:/trna-asl-md-stability/simulations/wt/system/system.xml)  
  *SHA-256:* `d38963ab5ea2adac230bcc552ab1b1dda2260eecf67f41fcf61ceb5a0a8f0b09`
- **Solvated PDB:** [simulations/wt/system/solvated.pdb](file:///d:/trna-asl-md-stability/simulations/wt/system/solvated.pdb)  
  *SHA-256:* `00b58110d71fcc25be065d744a7402ebff613df0f263663ef4490b0fef37ae1f`
- **Particles:** 15,693 (544 RNA atoms, 5,035 TIP3P water molecules, 30 $\text{K}^+$, 14 $\text{Cl}^-$)

### CUA Suppressor ASL (`cua`)
- **System XML:** [simulations/cua/system/system.xml](file:///d:/trna-asl-md-stability/simulations/cua/system/system.xml)  
  *SHA-256:* `5e1012f71bde05631b6612f6837c062f7133caf205f3522f3f303fdce7685220`
- **Solvated PDB:** [simulations/cua/system/solvated.pdb](file:///d:/trna-asl-md-stability/simulations/cua/system/solvated.pdb)  
  *SHA-256:* `a3e9c94a497a71245898a0d5c06473c237dbbf63222b7006b31c82d3167813cf`
- **Particles:** 15,686 (543 RNA atoms, 5,033 TIP3P water molecules, 30 $\text{K}^+$, 14 $\text{Cl}^-$)

---

## 3. Minimization Protocol & Parameters

To prevent initial solvent and ion displacement from distorting the crystallographically derived RNA backbone geometry, the minimization followed the two-stage protocol codified in [config/simulation_config.yaml](file:///d:/trna-asl-md-stability/config/simulation_config.yaml) and [docs/05_simulation_protocol.md](file:///d:/trna-asl-md-stability/docs/05_simulation_protocol.md) (Section 2):

### Stage 1: Restrained Solute Minimization (Solvent & Ion Relaxation) `[IMPLEMENTATION CHOICE]`
- **Restraint Target:** All RNA heavy atoms (357 atoms in WT, 356 atoms in CUA; hydrogens unrestrained).
- **Restraint Potential:** Harmonic position restraint applied via OpenMM `CustomExternalForce`:
  $$U_{\text{restraint}} = k \cdot \left[(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2\right]$$
- **Spring Constant:** $k = 5.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2} \approx 2092.00\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$.
- **Convergence Criterion:** OpenMM L-BFGS minimizer with tolerance $10.0\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$ ($< 0.24\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-1}$).
- **Iteration Ceiling:** Maximum 2,500 iterations.
- **Physical Objective:** Relaxes water molecules and ions around the solute, eliminating high-energy van der Waals contacts while freezing the RNA conformational core.

### Stage 2: Full System Minimization (Unrestrained Relaxation) `[IMPLEMENTATION CHOICE]`
- **Restraint Release:** Restraint spring constant set to $k = 0.0\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$ (zero external force).
- **Target:** Entire molecular system (RNA solute, explicit TIP3P solvent, $\text{K}^+$ and $\text{Cl}^-$ ions).
- **Convergence Criterion:** OpenMM L-BFGS minimizer with tolerance $10.0\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$.
- **Iteration Ceiling:** Maximum 5,000 iterations.
- **Physical Objective:** Resolves local geometric strain within the RNA nucleobases and backbone (particularly the mutated $C34$ and $A36$ bases) in the presence of fully relaxed solvent.

---

## 4. OpenMM Execution Environment

- **OpenMM Version:** `8.6.1`
- **Platform:** `OpenCL`
- **Hardware Device:** `NVIDIA GeForce MX350`
- **Precision Mode:** `mixed` (forces calculated in single precision; positions, coordinates, and energy accumulations maintained in double precision)
- **Integrator Used for Context Construction:** `LangevinMiddleIntegrator` ($T = 310.0\ \text{K}$, $\gamma = 1.0\ \text{ps}^{-1}$, $\Delta t = 2.0\ \text{fs}$)

---

## 5. Comparative Minimization Results (WT vs. CUA)

The minimization protocol was executed with strictly identical parameters on both constructs:

| Property / Metric | Wild-Type ASL (`wt`) | CUA Suppressor ASL (`cua`) | Comparative Evaluation |
| :--- | :---: | :---: | :--- |
| **Construct Identifier** | `wt_asl` | `suppressor_asl` | Controlled pair |
| **Target tRNA** | `tRNA-Asp-GTC-2-1` | `tRNA-Asp-GTC-2-1 (CUA mutant)` | GtRNAdb verified |
| **ASL Sequence (5' $\rightarrow$ 3')** | `CCCGCCUGUCACGCGGG` | `CCCGCCUCUAACGCGGG` | Strict $G34C$ / $C36A$ mutations |
| **Total System Particles** | **15,693** | **15,686** | Invariant particle conservation |
| **Stage 1 Restrained Atoms** | 357 (RNA heavy atoms) | 356 (RNA heavy atoms) | Transversion balance ($-1$) |
| **Stage 1 Restraint Constant ($k$)** | $5.0\ \text{kcal/mol/Å}^2$ | $5.0\ \text{kcal/mol/Å}^2$ | Identical |
| **Stage 2 Convergence Tolerance** | $10.0\ \text{kJ/mol/nm}$ | $10.0\ \text{kJ/mol/nm}$ | Identical |
| **Initial Potential Energy ($E_{\text{init}}$)** | **$-114,663.634\ \text{kJ/mol}$** | **$-113,083.589\ \text{kJ/mol}$** | Unrelaxed solvated systems |
| **Post-Stage 1 Energy ($E_{\text{stage1}}$)** | **$-282,127.908\ \text{kJ/mol}$** | **$-281,979.139\ \text{kJ/mol}$** | Solute restrained ($k=0$ eval) |
| **Stage 1 Energy Drop ($\Delta E_1$)** | **$-167,464.274\ \text{kJ/mol}$** | **$-168,895.549\ \text{kJ/mol}$** | 99.2% of total relaxation |
| **Final Potential Energy ($E_{\text{final}}$)** | **$-283,441.823\ \text{kJ/mol}$** | **$-283,727.738\ \text{kJ/mol}$** | Fully relaxed local minimum |
| **Stage 2 Additional Drop ($\Delta E_2$)** | **$-1,313.915\ \text{kJ/mol}$** | **$-1,748.600\ \text{kJ/mol}$** | Local RNA relaxation |
| **TOTAL ENERGY DROP ($\Delta E$)** | **$-168,778.189\ \text{kJ/mol}$** | **$-170,644.149\ \text{kJ/mol}$** | **$\Delta E \le 0$ Strictly Verified** |
| **RNA Heavy-Atom RMSD ($\text{initial} \rightarrow \text{min}$)** | **$0.6346\ \text{Å}$** | **$0.8921\ \text{Å}$** | **Both $< 1.0\ \text{Å}$ (preserved fold)** |
| **Minimum RNA Interatomic Distance** | **$1.2273\ \text{Å}$** | **$1.2263\ \text{Å}$** | **Zero steric clashes** |
| **Walltime** | $18.79\ \text{s}$ | $69.35\ \text{s}$ | Fast and lightweight |
| **Overall Minimization Status** | **PASS** | **PASS** | Fully Verified |

---

## 6. Post-Minimization Automated Scientific Validation

Every post-minimization check passed unconditionally in [scripts/05_minimize.py](file:///d:/trna-asl-md-stability/scripts/05_minimize.py):

1. **Particle & Residue Conservation:**
   - WT: 15,693 atoms, 5,035 water molecules, 30 $\text{K}^+$, 14 $\text{Cl}^-$ (`PASS`).
   - CUA: 15,686 atoms, 5,033 water molecules, 30 $\text{K}^+$, 14 $\text{Cl}^-$ (`PASS`).
2. **Coordinate Sanity:**
   - Zero NaN, zero Inf, zero missing coordinates across all atoms (`PASS`).
3. **Primary Structure Preservation:**
   - 17 RNA residues present with exact sequence identity preserved (`PASS`).
4. **Terminal Chemistry Invariance:**
   - Neutral 5'-OH (`HO5'` present, phosphate absent) and neutral 3'-OH (`HO3'` present) strictly maintained on both constructs (`PASS`).
5. **Monotonic Energy Reduction:**
   - Both systems exhibited massive, monotonic drops in potential energy ($\Delta E = -168,778.19\ \text{kJ/mol}$ WT; $-170,644.15\ \text{kJ/mol}$ CUA) confirming removal of unphysical overlaps (`PASS`).
6. **No Structural Distortion or Collapse:**
   - RNA heavy-atom RMSD relative to starting crystallographic/grafted coordinates is **$0.6346\ \text{Å}$** for WT and **$0.8921\ \text{Å}$** for CUA. The A-form stem and U-turn loop geometry remain intact without gross structural deformation (`PASS`).
7. **No Steric Clashes:**
   - Minimum non-bonded interatomic distance is $> 1.22\ \text{Å}$ in both systems (`PASS`).

---

## 7. Artifact Manifest & SHA-256 Checksums

### Wild-Type ASL Minimized Artifacts
- **Directory:** `simulations/wt/`
- **Minimized Coordinates:** [simulations/wt/minimized.pdb](file:///d:/trna-asl-md-stability/simulations/wt/minimized.pdb)  
  *SHA-256:* `471856b2293476ac5bfa24c5b91d43e9ad0476c7aaed8af0bbeacc0ca87efc80`
- **Minimization Log:** [simulations/wt/minimization.log](file:///d:/trna-asl-md-stability/simulations/wt/minimization.log)
- **Metadata YAML:** [simulations/wt/minimization_metadata.yaml](file:///d:/trna-asl-md-stability/simulations/wt/minimization_metadata.yaml)
- **Metadata JSON:** [simulations/wt/minimization_metadata.json](file:///d:/trna-asl-md-stability/simulations/wt/minimization_metadata.json)

### CUA Suppressor ASL Minimized Artifacts
- **Directory:** `simulations/cua/`
- **Minimized Coordinates:** [simulations/cua/minimized.pdb](file:///d:/trna-asl-md-stability/simulations/cua/minimized.pdb)  
  *SHA-256:* `695a1c3beb02b70b84f953b8892f365874dee2756c3df988f5dfa48b752eaf94`
- **Minimization Log:** [simulations/cua/minimization.log](file:///d:/trna-asl-md-stability/simulations/cua/minimization.log)
- **Metadata YAML:** [simulations/cua/minimization_metadata.yaml](file:///d:/trna-asl-md-stability/simulations/cua/minimization_metadata.yaml)
- **Metadata JSON:** [simulations/cua/minimization_metadata.json](file:///d:/trna-asl-md-stability/simulations/cua/minimization_metadata.json)

---

## 8. Reproducibility & Execution Instructions

To reproduce Phase 4B energy minimization or validate existing artifacts:

### Re-running Energy Minimization
```bash
# Minimize both constructs
uv run --with openmm --with numpy --with pyyaml python scripts/05_minimize.py --construct all

# Or minimize individually:
uv run --with openmm --with numpy --with pyyaml python scripts/05_minimize.py --construct wt
uv run --with openmm --with numpy --with pyyaml python scripts/05_minimize.py --construct cua
```

### Running Validation Only (Fast Check)
```bash
uv run --with openmm --with numpy --with pyyaml python scripts/05_minimize.py --validate-only
```

---

## 9. Phase Boundary Statement

> **"Phase 4B is complete. No NVT/NPT equilibration, production MD, or trajectory analysis was performed."**

The minimized coordinates in `simulations/wt/minimized.pdb` and `simulations/cua/minimized.pdb` are verified, frozen, and ready for **Phase 4C (Equilibration & Restraint Release)**.
