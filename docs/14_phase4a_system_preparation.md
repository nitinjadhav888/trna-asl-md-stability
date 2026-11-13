# 14 - Phase 4A Solvation & MD System Preparation Report

**Study Title:** tRNA ASL Structural Stability Study  
**Phase:** Phase 4A (Solvation & MD System Preparation)  
**Status:** Complete — Final Verdict: **PASS**  
**Execution Date:** November 13, 2026  
**Pipeline Script:** [scripts/04_prepare_system.py](file:///d:/trna-asl-md-stability/scripts/04_prepare_system.py)  

---

## 1. Scientific Purpose & Scope

The objective of Phase 4A is to construct fully solvated, electrically neutral, physiologically ion-equilibrated periodic simulation systems for both the Wild-Type 17-nt ASL and the CUA Suppressor 17-nt ASL, using validated molecular mechanics force fields prior to energy minimization and molecular dynamics.

### Strict Phase Boundary Enforcement
In strict compliance with the project roadmap:
- **Phase 4A Scope:** System assembly, explicit solvation, counterion neutralization, salt addition, force field parameterization, geometry audit, and serialization.
- **Explicit Exclusions:** **Zero energy minimization (Phase 4B), zero NVT/NPT equilibration (Phase 4C), zero production MD (Phase 5), and zero trajectory analysis** were performed. No conclusions regarding dynamic stability or biological suppression efficiency are drawn in this phase.

---

## 2. Molecular Modeling & Parameterization Protocol

### 2.1 Force Field Framework
- **RNA Solute:** AMBER OL3 (`ff99bsc0chiOL3`, encoded in OpenMM via `amber14/RNA.OL3.xml`).
  - Incorporates the Cornell et al. (1995) / Cheatham et al. (1999) `ff99` base parameters with the Pérez et al. (2007) `bsc0` $\alpha/\gamma$ backbone correction and the Banáš et al. (2010) `chiOL3` glycosidic torsion correction.
- **Explicit Water:** TIP3P (`amber14/tip3p.xml`).
  - Standard 3-site rigid water model with planar geometry ($r_{\text{OH}} = 0.9572\ \text{Å}$, $\theta_{\text{HOH}} = 104.52^\circ$, $q_{\text{O}} = -0.834\ e$, $q_{\text{H}} = +0.417\ e$).
- **Monovalent Ions:** Joung & Cheatham (2008) parameters for TIP3P (`ionsjc_tip3p`, encoded in `amber14/tip3p.xml`):
  - Potassium ($\text{K}^+$): $\sigma = 3.03796\ \text{Å}$ ($0.303796\ \text{nm}$), $\epsilon = 0.19368\ \text{kcal/mol}$ ($0.810369\ \text{kJ/mol}$).
  - Chloride ($\text{Cl}^-$): $\sigma = 4.47766\ \text{Å}$ ($0.447766\ \text{nm}$), $\epsilon = 0.03559\ \text{kcal/mol}$ ($0.148913\ \text{kJ/mol}$).
  - Strictly coupled to the TIP3P solvent model to prevent spurious ion clustering or artificial crystal nucleation.

### 2.2 Terminal Chemistry Verification
In accordance with SDR-09:
- **5'-Terminus (C27):** Neutral hydroxyl group ($5'\text{-OH}$), parameterized using OpenMM terminal template `C5`. Hydroxyl hydrogen `HO5'` is present; 5'-phosphate group ($P, OP1, OP2$) is absent.
- **3'-Terminus (G43):** Neutral hydroxyl group ($3'\text{-OH}$), parameterized using OpenMM terminal template `G3`. Terminal hydrogen `HO3'` is present.
- **Polyanionic Solute Net Charge:** 17 ribonucleotides contain exactly 16 internal phosphodiester links (each bearing a $-1.0000\ e$ charge). The unsolvated solute net charge is **$-16.0000\ e$** for both WT and CUA.

---

## 3. Box Geometry & Solute Clearance Verification

### 3.1 Selection of Cubic Periodic Boundary
An orthogonal **cubic periodic box** ($L = 5.500\ \text{nm} = 55.00\ \text{Å}$) was selected based on the following scientific criteria:
1. **Isotropic Rotational Clearance:** The isolated 17-nt ASL hairpin is an elongated molecule with a maximum dimension of $\approx 3.26\ \text{nm}$. A cubic box guarantees that even if the hairpin rotates or tumbles in solvent during unrestrained MD, the minimum clearance to any periodic face remains strictly $\ge 1.0\ \text{nm}$ ($10.0\ \text{Å}$) in all three Cartesian dimensions.
2. **Analysis Tool Compatibility:** Cubic boxes have strictly diagonal lattice vectors ($\mathbf{a} = (L,0,0), \mathbf{b} = (0,L,0), \mathbf{c} = (0,0,L)$), eliminating non-orthogonal coordinate transformation ambiguities across BioPython, MDAnalysis, OpenMM, and visualization engines.
3. **Identical Bounding Volumes:** Using an identical $55.00\ \text{Å}$ cubic box for both WT and CUA eliminates box volume discrepancies as a confounding variable in comparative analysis.

### 3.2 Programmatic Solute-to-Edge Clearance Audit
The clearance from every solute RNA atom $i$ to the nearest box boundary plane along each dimension $d \in \{x, y, z\}$ was computed as:
$$d_{i,d} = \min\left(r_{i,d} - (c_d - L/2),\ (c_d + L/2) - r_{i,d}\right)$$
where $c_d = \frac{1}{2}\left(\min_j r_{j,d} + \max_j r_{j,d}\right)$ is the geometric center of the solute.

- **Wild-Type ASL:**
  - Solute extent: $\Delta x = 2.101\ \text{nm}$, $\Delta y = 3.258\ \text{nm}$, $\Delta z = 2.109\ \text{nm}$
  - **Minimum solute clearance:** **$1.121\ \text{nm}$ ($11.21\ \text{Å}$)** $\ge 1.0\ \text{nm}$ `[PASS]`
- **CUA Suppressor ASL:**
  - Solute extent: $\Delta x = 2.101\ \text{nm}$, $\Delta y = 3.244\ \text{nm}$, $\Delta z = 2.109\ \text{nm}$
  - **Minimum solute clearance:** **$1.128\ \text{nm}$ ($11.28\ \text{Å}$)** $\ge 1.0\ \text{nm}$ `[PASS]`

---

## 4. Neutralization & Salt Concentration Protocol

### 4.1 Strict Conceptual Separation: Counterions vs. Added Salt
In accordance with Section 4 of the project brief, ions are explicitly divided into two independent chemical classes:

1. **Neutralizing Counterions ($16\ \text{K}^+$):**
   - The 16 internal phosphodiester groups carry a combined charge of $-16.0000\ e$.
   - Exactly $16\ \text{K}^+$ counterions are placed to neutralize the RNA backbone, bringing the polyanion net charge to zero.
2. **Physiological Added Salt ($14\ \text{K}^+$ and $14\ \text{Cl}^-$):**
   - Added to simulate a physiological ionic strength of approximately $0.15\ \text{M}$ ($150\ \text{mM}$) KCl.
   - Consists of exactly 14 stoichiometric neutral ion pairs ($14\ \text{K}^+$ cations and $14\ \text{Cl}^-$ anions).
3. **Total Ion Inventory:**
   - Total $\text{K}^+$ cations: $16 + 14 = \mathbf{30\ \text{K}^+}$
   - Total $\text{Cl}^-$ anions: $\mathbf{14\ \text{Cl}^-}$
   - Total System Net Charge: $30(+1) + 14(-1) - 16 = \mathbf{0.0000\ e}$ `[PASS]`

### 4.2 Exact Salt Concentration Analysis
Salt concentration was evaluated under both standard physical definitions:

1. **Box-Volume Concentration ($C_{\text{box}}$):**
   $$C_{\text{box}} = \frac{N_{\text{pairs}}}{N_A \times V_{\text{box}}} = \frac{14}{(6.02214 \times 10^{23}\ \text{mol}^{-1}) \times (1.66375 \times 10^{-22}\ \text{L})} = \mathbf{0.1397\ \text{M}}\ (139.7\ \text{mM})$$
   - Deviation from $0.150\ \text{M}$ target: **$-6.85\%$**
2. **Solvent-Volume Concentration ($C_{\text{solvent}}$):**
   Based on the standard partial molar volume of pure water at 310 K ($55.4\ \text{mol/L}$):
   $$V_{\text{water}} = \frac{N_{\text{water}}}{N_A \times 55.4\ \text{mol/L}} \approx 1.509 \times 10^{-22}\ \text{L}$$
   $$C_{\text{solvent}} = \frac{14}{(6.02214 \times 10^{23}\ \text{mol}^{-1}) \times (1.509 \times 10^{-22}\ \text{L})} = \mathbf{0.1540\ \text{M}}\ (154.0\ \text{mM})$$
   - Deviation from $0.150\ \text{M}$ target: **$+2.69\%$**

Both metrics demonstrate that the added salt rigorously brackets the $150\ \text{mM}$ target within physiological experimental tolerance ($<7\%$).

---

## 5. Comparative System Summary (WT vs. CUA)

| System Property | Wild-Type ASL (`wt_asl`) | CUA Suppressor ASL (`suppressor_asl`) | Comparative Assessment |
| :--- | :---: | :---: | :--- |
| **Primary Sequence (5' $\rightarrow$ 3')** | `CCCGCCUGUCACGCGGG` | `CCCGCCUCUAACGCGGG` | Strict $G34C$ / $C36A$ mutations |
| **Sequence Length** | 17 ribonucleotides | 17 ribonucleotides | Identical |
| **Sprinzl Coordinate Range** | Residues 27 to 43 | Residues 27 to 43 | Identical |
| **Terminal Chemistry** | Neutral 5'-OH, 3'-OH | Neutral 5'-OH, 3'-OH | Identical |
| **RNA Heavy Atoms** | 357 | 356 | Transversion balance ($-1$) |
| **RNA Total Atoms (with H)** | 544 | 543 | Transversion balance ($-1$) |
| **Unsolvated RNA Net Charge** | $-16.0000\ e$ | $-16.0000\ e$ | Identical |
| **Water Model** | TIP3P explicit solvent | TIP3P explicit solvent | Identical |
| **Water Molecules ($N_{\text{H2O}}$)** | 5,035 ($15,105$ atoms) | 5,033 ($15,099$ atoms) | Difference of only 2 waters (0.04%) |
| **Neutralizing Counterions** | $16\ \text{K}^+$ | $16\ \text{K}^+$ | Identical |
| **Added Salt Pairs (KCl)** | $14\ \text{K}^+ / 14\ \text{Cl}^-$ | $14\ \text{K}^+ / 14\ \text{Cl}^-$ | Identical |
| **Total $\text{K}^+$ Ions** | 30 | 30 | Identical |
| **Total $\text{Cl}^-$ Ions** | 14 | 14 | Identical |
| **Total System Particles** | **15,693** | **15,686** | Difference of 7 atoms (0.04%) |
| **Box Geometry** | Cubic (orthogonal) | Cubic (orthogonal) | Identical |
| **Box Length ($L$)** | $5.5000\ \text{nm}$ ($55.00\ \text{Å}$) | $5.5000\ \text{nm}$ ($55.00\ \text{Å}$) | Identical |
| **Box Volume ($V_{\text{box}}$)** | $166.375\ \text{nm}^3$ | $166.375\ \text{nm}^3$ | Identical |
| **Minimum Solute Clearance** | **$1.121\ \text{nm}$ ($11.21\ \text{Å}$)** | **$1.128\ \text{nm}$ ($11.28\ \text{Å}$)** | Both $> 1.0\ \text{nm}$ |
| **Salt Molarity (Solvent Vol)** | $0.1540\ \text{M}$ ($154.0\ \text{mM}$) | $0.1541\ \text{M}$ ($154.1\ \text{mM}$) | Difference $< 0.1\ \text{mM}$ |
| **Salt Molarity (Box Vol)** | $0.1397\ \text{M}$ ($139.7\ \text{mM}$) | $0.1397\ \text{M}$ ($139.7\ \text{mM}$) | Identical |
| **Final System Net Charge** | **$0.000000\ e$** | **$0.000000\ e$** | Exactly neutral |
| **OpenMM System Forces** | 5 (Bond, Angle, Torsion, Nonbonded, CMMotion) | 5 (Bond, Angle, Torsion, Nonbonded, CMMotion) | Identical |
| **Preparation Status** | **PASS** | **PASS** | Fully Verified |

---

## 6. Verification Checklist & Validation Results

Every system integrity requirement was programmatically checked during pipeline execution:

- [x] **Structural Identity:** Exactly 17 RNA residues present in both constructs (`PASS`).
- [x] **Sequence Invariance:** WT sequence matches `CCCGCCUGUCACGCGGG`; CUA matches `CCCGCCUCUAACGCGGG` (`PASS`).
- [x] **Terminal Chemistry Preserved:** Neutral 5'-OH (`HO5'` present, no phosphate) and neutral 3'-OH (`HO3'` present) on both (`PASS`).
- [x] **Force Field Parameterization:** AMBER OL3 + TIP3P loaded with zero missing templates or unparameterized atoms (`PASS`).
- [x] **Solute Clearance:** Programmatically verified minimum distance to box boundary is $> 1.0\ \text{nm}$ ($11.21\ \text{Å}$ WT, $11.28\ \text{Å}$ CUA) (`PASS`).
- [x] **Unsolvated Net Charge:** Exactly $-16.0000\ e$ (`PASS`).
- [x] **Counterion Neutralization:** Exactly $16\ \text{K}^+$ counterions placed to neutralize the polyanion (`PASS`).
- [x] **Salt Addition:** Exactly 14 neutral KCl pairs placed, achieving $0.154\ \text{M}$ solvent-volume concentration (`PASS`).
- [x] **Final System Net Charge:** Exactly $0.000000\ e$ for both systems (`PASS`).
- [x] **Coordinate Sanity:** Zero NaN, zero Inf, zero duplicate atoms, zero unphysical coordinates (`PASS`).
- [x] **Determinism:** Seed `42` enforced for exact algorithmic replay (`PASS`).

---

## 7. Artifact Manifest & SHA-256 Checksums

### Wild-Type ASL System Artifacts
- **Directory:** `simulations/wt/system/`
- **Solvated Structure:** [simulations/wt/system/solvated.pdb](file:///d:/trna-asl-md-stability/simulations/wt/system/solvated.pdb)  
  *SHA-256:* `00b58110d71fcc25be065d744a7402ebff613df0f263663ef4490b0fef37ae1f`
- **OpenMM System XML:** [simulations/wt/system/system.xml](file:///d:/trna-asl-md-stability/simulations/wt/system/system.xml)  
  *SHA-256:* `d38963ab5ea2adac230bcc552ab1b1dda2260eecf67f41fcf61ceb5a0a8f0b09`
- **Metadata JSON:** [simulations/wt/system/system_metadata.json](file:///d:/trna-asl-md-stability/simulations/wt/system/system_metadata.json)
- **Preparation Log:** [simulations/wt/system/preparation.log](file:///d:/trna-asl-md-stability/simulations/wt/system/preparation.log)

### CUA Suppressor ASL System Artifacts
- **Directory:** `simulations/cua/system/`
- **Solvated Structure:** [simulations/cua/system/solvated.pdb](file:///d:/trna-asl-md-stability/simulations/cua/system/solvated.pdb)  
  *SHA-256:* `a3e9c94a497a71245898a0d5c06473c237dbbf63222b7006b31c82d3167813cf`
- **OpenMM System XML:** [simulations/cua/system/system.xml](file:///d:/trna-asl-md-stability/simulations/cua/system/system.xml)  
  *SHA-256:* `5e1012f71bde05631b6612f6837c062f7133caf205f3522f3f303fdce7685220`
- **Metadata JSON:** [simulations/cua/system/system_metadata.json](file:///d:/trna-asl-md-stability/simulations/cua/system/system_metadata.json)
- **Preparation Log:** [simulations/cua/system/preparation.log](file:///d:/trna-asl-md-stability/simulations/cua/system/preparation.log)

---

## 8. Known Limitations & Dependencies for Phase 4B

### Inherent System Limitations
1. **Unrelaxed Starting Solvent-Solute Contacts:** Water molecules and ions were placed around rigid solute coordinates. Minor steric overlap or sub-optimal hydrogen-bonding orientations between solvent and RNA are present and require multi-stage energy minimization in Phase 4B.
2. **Unrelaxed Mutated Nucleobases:** The grafted $C34$ and $A36$ bases have ideal rigid-body bond lengths and angles, but their local van der Waals packing and electrostatic interactions with adjacent loop residues have not been energy-minimized.
3. **Rigid Water Model:** TIP3P is a non-polarizable model; charge polarization effects from the polyanionic RNA backbone are represented implicitly through fixed partial charges.

### Dependencies for Phase 4B (Energy Minimization)
Phase 4B requires:
1. Loading the serialized OpenMM systems (`system.xml`) and coordinates (`solvated.pdb`).
2. **Stage 1 (Restrained Solute Minimization):** Harmonic position restraints ($k = 5.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2} \approx 2092\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$) applied to all RNA heavy atoms while water and ions relax (2,500 steps steepest descent + 2,500 steps L-BFGS).
3. **Stage 2 (Full System Unrestrained Minimization):** Restraints released; entire system minimized until maximum force $F_{\text{max}} < 10.0\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$.
4. Output minimized coordinates (`minimized.pdb`) and state XML.

---

## 9. Final Phase 4A Status

> **"Phase 4A is complete. No minimization, equilibration, or production MD was performed."**

The repository is fully verified and ready for **Phase 4B only**.

