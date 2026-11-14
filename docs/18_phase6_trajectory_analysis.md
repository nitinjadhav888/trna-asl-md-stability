# Phase 6: Trajectory Analysis and Structural Comparison

## Executive Summary & Phase Boundary Scope
This document records the complete trajectory analysis, structural observables, replicate-level statistics, and sampling convergence assessments for **Phase 6** of the `trna-asl-md-stability` project.

> [!IMPORTANT]
> **Strict Phase Boundary Adherence:**
> This phase performs **trajectory-based structural analysis only**.
> - It does **NOT** perform biological interpretation of cellular suppression efficiency.
> - It does **NOT** claim that an MD structural difference proves suppression or reporter activity.
> - It does **NOT** write the final $\le 4$-page report (reserved for Phase 7).
> - All six production trajectories ($120.0\ \text{ns}$ aggregate sampling across $N=3$ WT and $N=3$ CUA independent simulation replicates) remained strictly immutable and unmodified.

---

## 1. Data Integrity and Trajectory Audit
Prior to structural analysis, an automated repository audit was conducted via `scripts/08_audit_production_data.py`. All six production directories were audited:

| Replicate Tag | Construct | Replicate Index | PRNG Seed | Trajectory File | Net Steps | Simulated Duration | Frame Count | Solute Atoms | QC Audit Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `wt_rep1` | WT | 1 | 5101 | `trajectory.dcd` | 10,000,000 | 20.00 ns | 2,000 | 544 | **PASS** (100%) |
| `wt_rep2` | WT | 2 | 5102 | `trajectory.dcd` | 10,000,000 | 20.00 ns | 2,000 | 544 | **PASS** (100%) |
| `wt_rep3` | WT | 3 | 5103 | `trajectory.dcd` | 10,000,000 | 20.00 ns | 2,000 | 544 | **PASS** (100%) |
| `cua_rep1` | CUA | 1 | 5201 | `trajectory.dcd` | 10,000,000 | 20.00 ns | 2,000 | 543 | **PASS** (100%) |
| `cua_rep2` | CUA | 2 | 5202 | `trajectory.dcd` | 10,000,000 | 20.00 ns | 2,000 | 543 | **PASS** (100%) |
| `cua_rep3` | CUA | 3 | 5203 | `trajectory.dcd` | 10,000,000 | 20.00 ns | 2,000 | 543 | **PASS** (100%) |

All raw production trajectories and checkpoints remain safely excluded from Git tracking via `.gitignore`.

---

## 2. Analysis Data Preparation & Workflow
To ensure fully reproducible and artifact-free trajectory coordinates:
- **Periodic Boundary Condition (PBC) Imaging:** Implemented in `scripts/09_extract_trajectory.py` using MDTraj `image_molecules()` to eliminate periodic boundary jumping and wrap water/ions while preserving the intact covalent topology of the 17-nt solute.
- **Analysis-Ready Trajectories:** Written to `simulations/{replicate_tag}/trajectory/solute_imaged.dcd` (~12.6 MB per replicate; 2,000 frames each; saving stride 10.0 ps).
- **Reference Coordinates:** Unified solute reference structures extracted from Phase 4C final equilibrated structures (`solute_ref.pdb`).
- **Software Stack:** MDTraj 1.10+, NumPy 2.2+, Pandas 2.2+, SciPy 1.15+, Matplotlib 3.10+, Seaborn 0.13+.
- **Frame Selection:** No initial equilibration frames were discarded arbitrarily. The entire 20.0 ns (2,000 frames) was analyzed uniformly for all six replicates; final-window properties (10.0–20.0 ns, frames 1,000–2,000) are explicitly tabulated alongside full-trajectory statistics.

---

## 3. Structural Metric Definitions & Specifications

### 3.1 Backbone Root-Mean-Square Deviation (RMSD)
Backbone atom selection: Six core atoms per nucleotide: `P`, `O5'`, `C5'`, `C4'`, `C3'`, `O3'`.
- **Full ASL RMSD:** Trajectory superposed on reference using all 102 backbone atoms (17 residues $\times$ 6 atoms); RMSD evaluated across all 102 atoms.
- **5-bp Stem RMSD:** Trajectory superposed on reference using stem backbone atoms (residues 27–31 and 39–43, 60 atoms); RMSD evaluated across these 60 atoms.
- **7-nt Loop Global RMSD (Stem-Aligned):** Trajectory remains superposed on the rigid 5-bp stem backbone reference; Cartesian RMS displacement evaluated on the 7 loop residues (positions 32–38, 42 atoms) to capture whole-loop rigid-body reorientation, tilt, and domain hinge motions.
- **7-nt Loop Internal RMSD (Loop-Aligned):** Trajectory superposed directly on the loop backbone atoms (positions 32–38); RMSD evaluated on these 42 atoms to isolate internal loop conformational rearrangements from stem hinge tilting.

### 3.2 Per-Residue Root-Mean-Square Fluctuation (RMSF)
Trajectories were aligned to the 5-bp stem backbone reference frame. Fluctuation around the average structure was computed for all atoms using `md.rmsf()`. Per-residue backbone RMSF was calculated as the root-mean-square of the atomic fluctuations of the six backbone atoms for each residue.

### 3.3 Stem Base-Pair Hydrogen-Bond Occupancy
Canonical Watson-Crick G-C base pairs:
1. `C27:G43` (5'-terminal pair)
2. `C28:G42`
3. `C29:G41`
4. `G30:C40`
5. `C31:G39` (loop-closing pair)

Geometric criteria for each hydrogen bond:
- Heavy-atom donor-to-acceptor distance: $d(\text{Donor}, \text{Acceptor}) \le 3.5\ \text{Å}$.
- Hydrogen-to-acceptor distance: $d(\text{H}, \text{Acceptor}) \le 2.5\ \text{Å}$.
- Angle: $\theta(\text{Donor}-\text{H}\cdots\text{Acceptor}) \ge 120.0^\circ$.
- Specific G-C hydrogen bonds evaluated:
  - **HB1:** G(N1)-H1 $\cdots$ C(N3)
  - **HB2:** G(N2)-H21/H22 $\cdots$ C(O2)
  - **HB3:** C(N4)-H41/H42 $\cdots$ G(O6)
- **Base Pair Intact Criterion:** $\ge 2$ of the 3 canonical Watson-Crick hydrogen bonds formed simultaneously.

### 3.4 Canonical Loop U-Turn Diagnostic Contact
Evaluates the canonical tertiary interaction that stabilizes the anticodon loop fold:
- **Donor:** U33(N3) (with bound proton H3).
- **Target Phosphate:** Non-bridging phosphate oxygens (`OP1`, `OP2`) of position 36:
  - WT: C36 phosphate (`resid` 9, C36 in WT).
  - CUA: A36 phosphate (`resid` 9, A36 suppressor mutation).
- **Contact Metric:** Minimum heavy-atom distance $d_{\text{min}} = \min(d(\text{N3}, \text{OP1}), d(\text{N3}, \text{OP2}))$.
- **Intact H-Bond Criterion:** $d_{\text{min}} \le 3.5\ \text{Å}$, $d(\text{H3}, \text{OP}) \le 2.5\ \text{Å}$, and $\angle(\text{N3}-\text{H3}\cdots\text{OP}) \ge 120.0^\circ$.

### 3.5 3'-Anticodon Base Stacking Geometry
Evaluated across adjacent nucleobases in the 3'-anticodon stack:
- Step 34–35 (Anticodon base 1 to base 2; G34-U35 in WT, C34-U35 in CUA)
- Step 35–36 (Anticodon base 2 to base 3; U35-C36 in WT, U35-A36 in CUA)
- Step 36–37 (Anticodon base 3 to invariant purine 37; C36-A37 in WT, A36-A37 in CUA)
- Step 37–38 (Purine 37 to closing nucleotide C38; A37-C38)

Geometric criteria (Singular Value Decomposition best-fit base planes and ring centroids):
- Pyrimidine ring atoms: `N1`, `C2`, `N3`, `C4`, `C5`, `C6`.
- Purine ring atoms: `N1`, `C2`, `N3`, `C4`, `C5`, `C6`, `N7`, `C8`, `N9`.
- **Standard Stacking:** Centroid distance $R_{\text{centroid}} \le 4.5\ \text{Å}$ AND inter-plane normal angle $\theta \le 30.0^\circ$.
- **Permissive Stacking:** $R_{\text{centroid}} \le 5.0\ \text{Å}$ AND $\theta \le 35.0^\circ$.

### 3.6 Additional Structural Metric: Nucleotide 34 Ribose Sugar Pucker
Endocyclic ribose furanose ring torsions:
- $\nu_0$: `C4'-O4'-C1'-C2'`
- $\nu_1$: `O4'-C1'-C2'-C3'`
- $\nu_2$: `C1'-C2'-C3'-C4'`
- $\nu_3$: `C2'-C3'-C4'-O4'`
- $\nu_4$: `C3'-C4'-O4'-C1'`

Pseudorotation phase angle $P$ (Altona & Sundaralingam, JACS 1972):
$$\tan(P) = \frac{(\nu_4 + \nu_1) - (\nu_3 + \nu_0)}{2 \cdot \nu_2 \cdot (\sin 36^\circ + \sin 72^\circ)},\quad P \in [0^\circ, 360^\circ)$$
- **North ($C3'\text{-endo}$, canonical A-form):** $0^\circ \le P < 36^\circ$ or $324^\circ \le P < 360^\circ$.
- **East ($O4'\text{-endo}$):** $36^\circ \le P < 144^\circ$.
- **South ($C2'\text{-endo}$, B-form-like):** $144^\circ \le P < 216^\circ$.
- **West ($C1'\text{-exo}/C4'\text{-endo}$):** $216^\circ \le P < 324^\circ$.

---

## 4. Replicate-Level Results

### 4.1 Master Replicate Data Table
Values reported for each independent 20.0 ns simulation run ($N=3$ WT, $N=3$ CUA):

| Observable | WT Replicate 1 | WT Replicate 2 | WT Replicate 3 | CUA Replicate 1 | CUA Replicate 2 | CUA Replicate 3 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Full ASL RMSD (Mean $\pm$ SD, Å)** | $1.845 \pm 0.336$ | $1.813 \pm 0.303$ | $2.139 \pm 0.453$ | $2.568 \pm 0.309$ | $2.293 \pm 0.502$ | $2.521 \pm 0.330$ |
| **Full ASL RMSD (Final 10 ns, Å)** | $1.909$ | $1.861$ | $2.187$ | $2.532$ | $2.589$ | $2.577$ |
| **Stem RMSD (Mean $\pm$ SD, Å)** | $0.870 \pm 0.199$ | $0.863 \pm 0.191$ | $0.990 \pm 0.242$ | $1.163 \pm 0.198$ | $1.121 \pm 0.211$ | $1.142 \pm 0.205$ |
| **Loop Global RMSD (Stem-Aligned, Å)** | $3.614 \pm 0.928$ | $3.421 \pm 0.852$ | $4.359 \pm 1.281$ | $5.569 \pm 1.075$ | $4.865 \pm 1.476$ | $5.324 \pm 1.064$ |
| **Loop Internal RMSD (Loop-Aligned, Å)**| $1.619 \pm 0.540$ | $1.799 \pm 0.381$ | $1.635 \pm 0.320$ | $2.250 \pm 0.267$ | $1.813 \pm 0.571$ | $2.141 \pm 0.359$ |
| **Stem 5-bp Intact Occupancy (%)** | $99.14\%$ | $99.98\%$ | $98.70\%$ | $99.27\%$ | $99.44\%$ | $99.50\%$ |
| **U-Turn HB Occupancy (%)** | $97.00\%$ | $98.30\%$ | $67.50\%$ | $0.00\%$ | $0.05\%$ | $0.05\%$ |
| **U-Turn Contact ($\le 3.5\ \text{Å}$) (%)** | $97.70\%$ | $99.50\%$ | $70.30\%$ | $0.00\%$ | $0.30\%$ | $0.05\%$ |
| **U-Turn Distance ($d_{\text{min}}$, Å)** | $2.90 \pm 0.32$ | $2.87 \pm 0.15$ | $3.42 \pm 0.87$ | $5.27 \pm 0.57$ | $6.20 \pm 1.42$ | $6.01 \pm 1.92$ |
| **Stacking Step 34–35 Occupancy (%)** | $93.10\%$ | $95.45\%$ | $84.20\%$ | $2.05\%$ | $40.90\%$ | $14.80\%$ |
| **Stacking Step 35–36 Occupancy (%)** | $51.85\%$ | $50.75\%$ | $40.25\%$ | $50.65\%$ | $55.60\%$ | $50.85\%$ |
| **Stacking Step 36–37 Occupancy (%)** | $50.55\%$ | $50.45\%$ | $64.40\%$ | $59.15\%$ | $69.45\%$ | $59.70\%$ |
| **Stacking Step 37–38 Occupancy (%)** | $77.75\%$ | $90.40\%$ | $87.25\%$ | $98.75\%$ | $76.00\%$ | $96.70\%$ |
| **Pos 34 North ($C3'\text{-endo}$) Occ (%)** | $86.65\%$ | $84.90\%$ | $86.95\%$ | $77.80\%$ | $84.20\%$ | $79.55\%$ |
| **Pos 34 East ($O4'\text{-endo}$) Occ (%)** | $4.85\%$ | $3.20\%$ | $6.60\%$ | $22.15\%$ | $15.65\%$ | $20.30\%$ |
| **Pos 34 Median Phase Angle $P$ ($^\circ$)** | $331.2^\circ$ | $333.2^\circ$ | $329.7^\circ$ | $26.1^\circ$ | $25.3^\circ$ | $25.9^\circ$ |

---

## 5. Statistical Comparison (WT vs CUA)

### 5.1 Replicate-Level Hypothesis Testing
Following the strict mandate of Section 8, the three WT trajectories and three CUA trajectories were treated as the $N=3$ independent statistical units. Trajectory frames were **never** pooled to inflate degrees of freedom.
- **Statistical Test:** Two-sample unpaired Welch's $t$-test (unequal variance, Satterthwaite degrees of freedom).
- **Effect Size:** Cohen's $d$ (pooled sample standard deviation).
- **Multiple Comparison Control:** Benjamini-Hochberg False Discovery Rate ($q$-value) computed across all 15 structural metrics.

| Observable Metric | WT Mean $\pm$ SD ($N=3$) | WT 95% CI | CUA Mean $\pm$ SD ($N=3$) | CUA 95% CI | Difference $\Delta$ | Cohen's $d$ | Welch's $t$ | Welch $df$ | $p$-value | FDR $q$-value |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Full ASL RMSD (Å)** | $1.932 \pm 0.180$ | $[1.486, 2.378]$ | $2.461 \pm 0.147$ | $[2.095, 2.827]$ | $+0.529$ | $+3.22$ | $3.94$ | $3.85$ | $0.0182$ | **$0.0390$** |
| **Full ASL RMSD (Final 10 ns, Å)** | $1.986 \pm 0.176$ | $[1.549, 2.422]$ | $2.566 \pm 0.030$ | $[2.491, 2.641]$ | $+0.581$ | $+4.61$ | $5.64$ | $2.12$ | $0.0264$ | **$0.0439$** |
| **5-bp Stem RMSD (Å)** | $0.908 \pm 0.071$ | $[0.731, 1.085]$ | $1.142 \pm 0.021$ | $[1.091, 1.194]$ | $+0.234$ | $+4.46$ | $5.46$ | $2.34$ | $0.0224$ | **$0.0421$** |
| **Loop Global RMSD (Å)** | $3.798 \pm 0.496$ | $[2.567, 5.029]$ | $5.253 \pm 0.358$ | $[4.365, 6.141]$ | $+1.455$ | $+3.37$ | $4.12$ | $3.64$ | $0.0177$ | **$0.0390$** |
| **Loop Internal RMSD (Å)** | $1.685 \pm 0.100$ | $[1.437, 1.933]$ | $2.068 \pm 0.227$ | $[1.503, 2.633]$ | $+0.383$ | $+2.18$ | $2.67$ | $2.74$ | $0.0830$ | $0.1166$ |
| **Stem Base-Pair Occupancy** | $0.993 \pm 0.006$ | $[0.977, 1.009]$ | $0.994 \pm 0.001$ | $[0.991, 0.997]$ | $+0.001$ | $+0.29$ | $0.35$ | $2.13$ | $0.7563$ | $0.7563$ |
| **U-Turn HB Occupancy** | $0.876 \pm 0.174$ | $[0.443, 1.309]$ | $0.0003 \pm 0.0003$ | $[-0.0004, 0.0011]$ | $-0.876$ | $-7.11$ | $-8.71$ | $2.00$ | $0.0129$ | **$0.0390$** |
| **U-Turn Contact ($\le 3.5\ \text{Å}$) Occ** | $0.892 \pm 0.164$ | $[0.485, 1.298]$ | $0.0012 \pm 0.0016$ | $[-0.0028, 0.0052]$ | $-0.891$ | $-7.70$ | $-9.43$ | $2.00$ | $0.0111$ | **$0.0390$** |
| **U-Turn Distance ($d_{\text{min}}$, Å)** | $3.060 \pm 0.309$ | $[2.292, 3.828]$ | $5.827 \pm 0.490$ | $[4.609, 7.045]$ | $+2.767$ | $+6.75$ | $8.27$ | $3.37$ | $0.0024$ | **$0.0354$** |
| **Stacking Step 34–35 Occupancy** | $0.909 \pm 0.059$ | $[0.762, 1.057]$ | $0.193 \pm 0.198$ | $[-0.299, 0.684]$ | $-0.717$ | $-4.90$ | $-6.00$ | $2.36$ | $0.0178$ | **$0.0390$** |
| **Stacking Step 35–36 Occupancy** | $0.476 \pm 0.064$ | $[0.317, 0.635]$ | $0.524 \pm 0.028$ | $[0.454, 0.593]$ | $+0.047$ | $+0.96$ | $1.18$ | $2.74$ | $0.3313$ | $0.3822$ |
| **Stacking Step 36–37 Occupancy** | $0.551 \pm 0.080$ | $[0.352, 0.751]$ | $0.628 \pm 0.058$ | $[0.484, 0.772]$ | $+0.076$ | $+1.09$ | $1.34$ | $3.64$ | $0.2590$ | $0.3238$ |
| **Stacking Step 37–38 Occupancy** | $0.851 \pm 0.066$ | $[0.688, 1.015]$ | $0.905 \pm 0.126$ | $[0.592, 1.217]$ | $+0.053$ | $+0.53$ | $0.65$ | $3.02$ | $0.5604$ | $0.6004$ |
| **Pos 34 Sugar North ($C3'\text{-endo}$)** | $0.862 \pm 0.011$ | $[0.834, 0.889]$ | $0.805 \pm 0.033$ | $[0.723, 0.887]$ | $-0.056$ | $-2.29$ | $-2.81$ | $2.44$ | $0.0855$ | $0.1166$ |
| **Pos 34 Sugar East ($O4'\text{-endo}$)** | $0.049 \pm 0.017$ | $[0.007, 0.091]$ | $0.194 \pm 0.033$ | $[0.110, 0.277]$ | $+0.145$ | $+5.45$ | $6.68$ | $2.97$ | $0.0071$ | **$0.0390$** |

---

## 6. Sampling Sufficiency & Convergence Assessment

### 6.1 Block Analysis (Four Contiguous 5-ns Time Windows)
Each 20.0 ns trajectory was partitioned into four 5.0 ns windows:
- **Block 1:** 0.0–5.0 ns (frames 0–500)
- **Block 2:** 5.0–10.0 ns (frames 500–1,000)
- **Block 3:** 10.0–15.0 ns (frames 1,000–1,500)
- **Block 4:** 15.0–20.0 ns (frames 1,500–2,000)

Detailed window metrics across all replicates:
- **Stem Stability:** Across all four blocks in all six trajectories, stem backbone RMSD remained strictly stable ($0.82$ to $1.19\ \text{Å}$) and base-pair intact occupancy remained $>98.5\%$.
- **WT Dynamics & Replicate Heterogeneity:**
  - `wt_rep1` and `wt_rep2` exhibited stationary behavior throughout the full 20 ns, maintaining intact U-turn contacts ($91.0\%$ to $99.4\%$ HB occupancy) and robust 34–35 stacking ($88.8\%$ to $99.2\%$).
  - In contrast, `wt_rep3` displayed conformational drift in the loop: while Block 1 showed canonical U-turn geometry (HB occ $94.6\%$, distance $2.93\ \text{Å}$), the loop underwent spontaneous structural breathing in Blocks 3 and 4, with U-turn HB occupancy dropping to $60.4\%$ and $30.2\%$, and distance increasing to $4.22\ \text{Å}$.
- **CUA Conformational Dynamics:**
  - In `cua_rep1`, the U-turn contact was disrupted from the start ($0.0\%$ occupancy), and 34–35 stacking dropped from $7.4\%$ in Block 1 to $<0.5\%$ in Blocks 2–4.
  - In `cua_rep2` and `cua_rep3`, the 34–35 dinucleotide step began with transient stacking in Block 1 ($89.2\%$ in rep 2; $58.0\%$ in rep 3), but unstacked dramatically after ~8–12 ns (dropping to $0.0\%$ in rep 2 Blocks 3–4; $<0.8\%$ in rep 3 Blocks 2–4).

### 6.2 Rigorous Distinction of Stability vs Convergence
In accordance with Section 9 requirements, four hierarchical levels of simulation behavior are distinguished:
1. **Numerical Stability:** **Confirmed.** All six trajectories maintained stable integration with 2.0 fs time steps, conserved total energy during microcanonical checks, and exhibited no particle clashes or box boundary instabilities.
2. **Local Structural Stability:** **Confirmed.** The 5-bp stem scaffold (C27:G43 through C31:G39) remained locked in an intact Watson-Crick A-form double helix across all 120 ns of aggregate sampling.
3. **Sampling Convergence:** **Not Achieved for the Loop.** Because individual replicates exhibited transitions occurring between 5 and 15 ns (such as the unstacking in `cua_rep2` and the loop breathing in `wt_rep3`), the loop conformational states have not achieved stationary equilibrium within a single 20 ns window.
4. **Global Conformational Equilibrium:** **Not Claimed.** An aggregate sampling of $3 \times 20\ \text{ns}$ is sufficient to demonstrate that CUA destabilizes the canonical U-turn contact and disrupts 34–35 stacking, but is insufficient to fully map the equilibrium free-energy landscape or relative populations of all accessible loop substates.

---

## 7. Generated Figures Index
All figures are saved in the [figures/](file:///d:/trna-asl-md-stability/figures) directory at 300 DPI publication resolution:

1. [figures/figure1_rmsd.png](file:///d:/trna-asl-md-stability/figures/figure1_rmsd.png):
   - Multi-panel plot: (A) Full ASL RMSD time series; (B) 5-bp stem RMSD time series; (C) 7-nt loop global RMSD (stem-aligned); (D) Replicate mean distributions with individual replicate dots.
2. [figures/figure2_rmsf.png](file:///d:/trna-asl-md-stability/figures/figure2_rmsf.png):
   - Per-residue backbone RMSF for all 17 nucleotides (tRNA positions 27–43), highlighting 5'-stem, 7-nt loop, anticodon, and 3'-stem domains.
3. [figures/figure3_stem_hbond.png](file:///d:/trna-asl-md-stability/figures/figure3_stem_hbond.png):
   - Watson-Crick intact occupancy across all five stem pairs (C27:G43, C28:G42, C29:G41, G30:C40, C31:G39) for WT and CUA.
4. [figures/figure4_uturn_diagnostic.png](file:///d:/trna-asl-md-stability/figures/figure4_uturn_diagnostic.png):
   - Loop U-turn tertiary interaction: (A) Distance time series; (B) Probability density distributions; (C) Replicate-level intact hydrogen-bond occupancy.
5. [figures/figure5_stacking.png](file:///d:/trna-asl-md-stability/figures/figure5_stacking.png):
   - Dinucleotide base stacking occupancy across steps 34–35, 35–36, 36–37, and 37–38.
6. [figures/figure6_sugar_pucker.png](file:///d:/trna-asl-md-stability/figures/figure6_sugar_pucker.png):
   - Nucleotide 34 ribose sugar pucker: (A) Pseudorotation phase angle $P$ distribution; (B) Conformational state occupancy breakdown (North, East, South, West).

---

## 8. Summary of Analysis Scripts & Machine-Readable Tables

### 8.1 Analysis Scripts
- `scripts/08_audit_production_data.py`: Pre-analysis file and integrity audit.
- `scripts/09_extract_trajectory.py`: PBC wrap imaging and solute trajectory extraction.
- `scripts/10_rmsd.py`: Multi-domain backbone RMSD analysis.
- `scripts/11_rmsf.py`: Per-residue backbone and heavy-atom RMSF analysis.
- `scripts/12_hbond_analysis.py`: Stem base-pair Watson-Crick hydrogen-bond analysis.
- `scripts/13_uturn_analysis.py`: Canonical U-turn diagnostic interaction analysis.
- `scripts/14_stacking_analysis.py`: 3'-anticodon base stacking geometry analysis.
- `scripts/15_additional_metric.py`: Nucleotide 34 ribose sugar pucker analysis.
- `scripts/16_statistics.py`: Master table collation, Welch's $t$-test, FDR, and block convergence.
- `scripts/17_generate_figures.py`: Publication figure rendering.

### 8.2 Analysis Tables
All tables are stored in [analysis/tables/](file:///d:/trna-asl-md-stability/analysis/tables):
- `rmsd_replicates.csv` & `rmsd_summary.csv` & `rmsd_time_series.csv`
- `rmsf_replicates.csv` & `rmsf_summary.csv`
- `stem_hbond_replicates.csv` & `stem_hbond_summary.csv`
- `uturn_replicates.csv` & `uturn_summary.csv` & `uturn_time_series.csv`
- `stacking_replicates.csv` & `stacking_summary.csv` & `stacking_time_series.csv`
- `sugar_pucker_replicates.csv` & `sugar_pucker_summary.csv` & `sugar_pucker_time_series.csv`
- `replicate_master_table.csv`
- `statistical_comparison.csv`
- `sampling_assessment.csv`

---

## 9. Next Steps
Phase 6 structural trajectory analysis is concluded.
All outputs are fully prepared for **Phase 7** (Final Interpretation, Synthesis, and $\le 4$-Page Publication Report).
