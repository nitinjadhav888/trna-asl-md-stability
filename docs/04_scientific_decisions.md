# 04 - Scientific Decision Log

This log formally records all scientific and technical decisions, their rationales, alternatives considered, and governance status.

Taxonomy:
- `[HARD REQUIREMENT]` — Explicitly dictated by the take-home prompt or immutable biological laws.
- `[SCIENTIFIC DECISION]` — Biophysical, structural, or methodological choice grounded in scientific literature.
- `[IMPLEMENTATION CHOICE]` — Operational software, numerical setting, or algorithmic heuristic.
- `[ENGINEERING CHECK]` — Quality-assurance threshold or pipeline sanity check.
- `[OPEN QUESTION]` — Unresolved decision deferred to a specific future pipeline phase.

---

## Record Index

| ID | Title | Status | Classification |
| :--- | :--- | :--- | :--- |
| **SDR-01** | ASL Structural Boundary (17-nt, 5-bp Stem, 7-nt Loop) | **APPROVED** | `[HARD REQUIREMENT]` |
| **SDR-02** | Force Field Parameterization: AMBER OL3 (`ff99bsc0chiOL3`) | **APPROVED** | `[HARD REQUIREMENT]` |
| **SDR-03** | Thermodynamic Ensemble ($NPT$ at 310 K, 1 bar) | **APPROVED** | `[HARD REQUIREMENT]` |
| **SDR-04** | Sampling Budget (3 replicates × 20 ns per construct = 120 ns) | **APPROVED** | `[HARD REQUIREMENT]` |
| **SDR-05** | Ionic Environment ($\approx 0.15\ \text{M}\ \text{KCl}$ in Explicit Water) | **APPROVED** | `[HARD REQUIREMENT]` |
| **SDR-06** | Solvent Model & Matched Ion Parameter Selection (TIP3P baseline) | **BASELINE ESTABLISHED** | `[IMPLEMENTATION CHOICE]` |
| **SDR-07** | Target Human Cytoplasmic tRNA Selection (`tRNA-Asp-GTC-2-1`) | **APPROVED (FROZEN)** | `[HARD REQUIREMENT]` / `[IMPLEMENTATION CHOICE]` |
| **SDR-08** | Selection & Justification of Metric 6 | **DEFERRED TO PHASE 5** | `[OPEN QUESTION]` |
| **SDR-09** | ASL Terminal Chemistry (Neutral 5'-OH / 3'-OH) | **APPROVED & IMPLEMENTED** | `[SCIENTIFIC DECISION]` |
| **SDR-10** | Replicate-Aware Statistical Inference & Sample Size Constraints | **METHODOLOGY DEFINED** | `[SCIENTIFIC DECISION]` |
| **SDR-11** | 3D Structural Modeling Template (PDB 6UGG) & Base Grafting | **APPROVED & IMPLEMENTED** | `[SCIENTIFIC DECISION]` / `[IMPLEMENTATION CHOICE]` |
| **SDR-12** | Pre-MD Structural & Force Field Validation Standards | **APPROVED & IMPLEMENTED** | `[ENGINEERING CHECK]` |
| **SDR-13** | Two-Stage Solvated Energy Minimization Protocol | **APPROVED & IMPLEMENTED** | `[IMPLEMENTATION CHOICE]` |
| **SDR-14** | Staged NVT/NPT Restrained Equilibration Protocol | **APPROVED & IMPLEMENTED** | `[IMPLEMENTATION CHOICE]` |

---

## Detailed Scientific Decision Records

### SDR-01: ASL Structural Boundary Definition
- **Status:** APPROVED `[HARD REQUIREMENT]`
- **Context:** An anticodon stem-loop (ASL) within a canonical tRNA spans residues 27–43 in Sprinzl numbering.
- **Decision:** Model an isolated 17-nucleotide hairpin comprising a 5-bp Watson-Crick stem (base pairs 27:43, 28:42, 29:41, 30:40, 31:39) and a 7-nucleotide anticodon loop (residues 32–38).
- **Rationale:** Minimal autonomous structural unit representing the anticodon domain mandated by the take-home requirements.

---

### SDR-02: Force Field Parameterization (AMBER OL3)
- **Status:** APPROVED `[HARD REQUIREMENT]`
- **Context:** Empirical force fields require accurate glycosidic torsions to avoid unphysical RNA laddering or loop collapse.
- **Decision:** Utilize the AMBER OL3 parameter set (`ff99bsc0chiOL3`, encoded in OpenMM as `amber14/RNA.OL3.xml`).
- **Rationale:** OL3 specifically corrects RNA $\chi$ glycosidic torsion biases and represents the established standard for RNA hairpins.

---

### SDR-03: Thermodynamic Ensemble (NPT at 310 K, 1 bar)
- **Status:** APPROVED `[HARD REQUIREMENT]`
- **Context:** Human cytoplasmic tRNAs function in vivo at physiological body temperature.
- **Decision:** Simulate in the isothermal-isobaric ensemble ($NPT$) at $T = 310.0\ \text{K}$ and $P = 1.0\ \text{bar}$.
- **Rationale:** Matches physiological human body temperature ($37^\circ\text{C}$) and standard pressure.

---

### SDR-04: Sampling Budget & Replicate Control
- **Status:** APPROVED `[HARD REQUIREMENT]`
- **Context:** RNA hairpin trajectories can become trapped in local kinetic minima.
- **Decision:** Execute 3 statistically independent replicates per construct (WT and CUA suppressor), each running for a minimum of 20 ns production, totaling 6 simulations and 120 ns aggregate sampling. Assign initial velocities using distinct pseudo-random seeds.
- **Rationale:** Multiple independent trajectories provide an assessment of replicate variability that cannot be captured by a single trajectory.

---

### SDR-05: Ionic Environment (~0.15 M KCl)
- **Status:** APPROVED `[HARD REQUIREMENT]`
- **Context:** RNA is a polyanion requiring monovalent counterions to screen phosphate repulsions.
- **Decision:** Solvate in explicit water, neutralize the 16 negative charges with $\text{K}^+$, and add $\text{K}^+$ and $\text{Cl}^-$ co-ions to reach an ionic strength of approximately $0.15\ \text{M}$.
- **Rationale:** Matches physiological human cytoplasmic potassium concentrations ($\sim 140\text{--}150\ \text{mM}\ \text{K}^+$).

---

### SDR-06: Solvent Model & Matched Ion Parameter Selection
- **Status:** BASELINE ESTABLISHED `[IMPLEMENTATION CHOICE]`
- **Context:** The take-home prompt explicitly permits TIP3P or an equivalent model. The water model must be compatible with AMBER OL3 and matched with appropriate ion parameters.
- **Decision:** 
  - Establish **TIP3P** as the primary baseline water model for this project because AMBER OL3 (`ff99bsc0chiOL3`) was parameterized and validated natively in TIP3P.
  - Pair TIP3P with matched monovalent ion parameters (e.g., Joung-Cheatham `ionsjc_tip3p` or standard AMBER ion libraries in `amber14/tip3p.xml`).
  - Alternative water models (such as OPC or TIP4P-EW) remain optional only if accompanied by explicit, documented scientific justification and strictly matched ion parameters (e.g., `ionsjc_opc`).
- **Rationale:** Enforces internal consistency across force field, solvent, and ion parameters, avoiding arbitrary or unverified solvent preferences.

---

### SDR-07: Target Human Cytoplasmic tRNA Selection (tRNA-Asp-GTC-2-1)
- **Status:** APPROVED (PHASE 2 COMPLETED & FROZEN) `[HARD REQUIREMENT]` / `[IMPLEMENTATION CHOICE]`
- **Context:** The study requires selecting exactly one wild-type human cytoplasmic tRNA from GtRNAdb whose native anticodon can be converted to the amber suppressor anticodon `CUA` by **exactly two nucleotide substitutions** ($\text{Hamming}(\text{Anticodon}_{\text{WT}}, \text{CUA}) = 2$).
- **Evidence & Screening Execution:**
  - Evaluated all 432 high-confidence *Homo sapiens* cytoplasmic tRNA genes in GtRNAdb (hg38).
  - Computed exact Hamming distance to `CUA` across all 48 unique anticodons; 199 genes (21 anticodons) satisfied $\text{Hamming} = 2$ `[HARD REQUIREMENT]`.
  - Filtered for canonical 17-nt ASL secondary structure (`>>>>>.......<<<<<`), invariant U33, anticodon at loop positions 3–5 (residues 34–36), purine at position 6 (residue 37), and pure Watson-Crick 5-bp stem. Exactly 138 genes across 30 unique ASL sequence clusters passed all filters `[HARD REQUIREMENT]`.
- **Decision:** Select **Human Cytoplasmic tRNA-Asp-GTC-2-1** (GtRNAdb ID: `tRNA-Asp-GTC-2-1`, tRNAscan-SE ID: `chr1.trna92`, genomic locus `chr1:161440825-161440896`).
  - Native WT Anticodon: `5'-GUC-3'` (DNA: `GTC`).
  - Target Suppressor Anticodon: `5'-CUA-3'` (DNA: `CTA`).
  - Hamming Distance: **2** (Substitutions: `G34C` and `C36A`; middle position `U35` is identical).
  - 17-nt ASL Sequence (WT): `5'-CCCGCCUGUCACGCGGG-3'`
  - 17-nt ASL Sequence (Suppressor): `5'-CCCGCCUCUAACGCGGG-3'`
- **Scientific Rationale & Distinctions:**
  1. *High-GC Stem Composition `[IMPLEMENTATION CHOICE]` (5/5 C:G Base Pairs):* The 5-bp stem consists entirely of Watson-Crick C:G pairs. High-GC stem composition was selected as an implementation/design criterion intended to reduce the likelihood of stem fraying during short ASL-only simulations (RSK-SCI-03). This is not claimed as a proven $\Delta G$ or guarantee of helical preservation; actual stem stability remains an empirical outcome to be directly measured via base-pair hydrogen-bond occupancy (Metric 3) in Phase 5 trajectory analysis.
  2. *Canonical U-Turn Architecture `[HARD REQUIREMENT / BIOLOGICAL FACT]`:* Features invariant U33 and invariant purine A37 flanking the anticodon triplet.
  3. *High Biological Representation & Genomic Copy Count `[IMPLEMENTATION CHOICE]`:* 13 total `tRNA-Asp-GTC` genes exist in the high-confidence cytosolic set. Of these, 11 share the identical 72-nt mature sequence of `tRNA-Asp-GTC-2-1` (`2-1` through `2-11`), and 12 share the identical 17-nt ASL sequence (`CCCGCCUGUCACGCGGG`), confirming it as the predominant human aspartate tRNA multi-gene family.
  4. *Zero Introns `[HARD REQUIREMENT]`:* Non-intronic transcript with verified mature sequence from GtRNAdb.

---

### SDR-08: Selection & Justification of Metric 6
- **Status:** DEFERRED TO PHASE 5 `[OPEN QUESTION]`
- **Context:** The take-home prompt requires one additional justified structural metric alongside RMSD, RMSF, H-bonds, U-turn diagnostic, and 3'-stacking.
- **Decision:** No metric is prematurely locked in Phase 1-3. Formal selection and scientific justification will occur during Phase 5 analysis from the following candidates:
  - *Candidate A:* **Anticodon Loop Radius of Gyration ($R_{g,\text{loop}}$)** — Quantifies global loop expansion or compaction.
  - *Candidate B:* **Anticodon Solvent-Accessible Surface Area (SASA)** — Quantifies nucleobase solvent exposure vs. lumen flipping.
  - *Candidate C:* **Stem-Loop Interhelical Angle / Axis Tilt** — Quantifies coaxial bending relative to the stem helix.

---

### SDR-09: ASL Terminal Chemistry
- **Status:** APPROVED & IMPLEMENTED `[SCIENTIFIC DECISION]`
- **Context:** Simulating an excised 17-nt hairpin requires defining terminal groups at residues 27 and 43.
- **Decision:** Maintain neutral **5'-OH** and **3'-OH** termini.
- **Rationale:** Prevents placing a repulsive $-2e$ charge on residue 27 directly opposite residue 43, minimizing artificial terminal base-pair fraying (RSK-SCI-03).
- **Implementation in Phase 3:** In `structures/wt/initial.pdb` and `structures/cua/initial.pdb`, residue 27 begins with $O5'$ (phosphate atoms $P, OP1, OP2$ omitted), matching OpenMM force field template `C5` (adds $HO5'$). Residue 43 terminates with $O3'$, matching template `G3` (adds $HO3'$).

---

### SDR-10: Replicate-Aware Statistical Inference & Sample Size Constraints
- **Status:** METHODOLOGY DEFINED `[SCIENTIFIC DECISION]`
- **Context:** Trajectory frames sampled at 10 ps intervals exhibit high temporal autocorrelation. Furthermore, with $N=3$ independent replicates per construct, statistical power for traditional null-hypothesis testing is severely constrained.
- **Mathematical Reality:**
  - A two-sided Mann-Whitney $U$ test with $N_1 = 3, N_2 = 3$ has $\binom{6}{3} = 20$ possible rank arrangements. The minimum attainable two-sided $p$-value is $2/20 = 0.10$.
  - Therefore, claiming that classical significance testing can achieve $p < 0.05$ with $N=3$ is mathematically impossible.
- **Decision:**
  1. **Primary Reporting:** Center reporting on replicate-level effect sizes (Hedges' $g$ / Cohen's $d$, difference in replicate means/medians) with 95% bootstrap confidence intervals ($10,000$ iterations).
  2. **Explicit Disclosure:** Transparently state sample-size and power limitations in the final report.
  3. **Pseudoreplication Guard:** Trajectory frames will not be treated as independent biological replicates. Autocorrelation times ($\tau_{\text{int}}$) and block averages will be computed to assess intra-trajectory stability, while inter-replicate comparisons remain the primary statistical unit.

---

### SDR-11: 3D Structural Modeling Template and Glycosidic Base Grafting
- **Status:** APPROVED & IMPLEMENTED `[SCIENTIFIC DECISION]` / `[IMPLEMENTATION CHOICE]`
- **Context:** Phase 3 requires constructing valid 3D coordinates for WT and CUA suppressor 17-nt ASLs under strictly controlled comparative conditions.
- **Decision:** Select experimental X-ray crystal structure of unmodified tRNA(Asp) (PDB: **6UGG**, Chain A, residues 27–43, 1.95 Å resolution) as the structural template, with deterministic rigid-body glycosidic base grafting for stem mutations ($U29C, A41G$) and suppressor mutations ($G34C, C36A$).
- **Rationale & Audit Confirmation:**
  1. The 7-nt loop of 6UGG (`CUGUCAC`, residues 32–38) is 100% sequence-identical to human cytoplasmic `tRNA-Asp-GTC-2-1`, preserving authentic crystallographic U-turn geometry ($U33\ \text{N3} \cdots 36\ \text{OP2} = 3.52\ \text{Å}$) and canonical $\text{anti}$ glycosidic torsions ($\chi$).
  2. Chain A was specifically verified over Chain B because Chain A adopts the canonical folded U-turn conformation ($3.52\ \text{Å}$), whereas Chain B exhibits an extended non-canonical loop ($13.36\ \text{Å}$) stabilized by crystal lattice contacts (Chan et al., *RNA* 2020).
  3. The 5-bp stem of 6UGG has 4 of 5 pairs already identical; the single compensatory mutation ($U29C, A41G$) preserves standard A-form Watson-Crick hydrogen bonding ($2.74 - 2.82\ \text{Å}$).
  4. Non-mutated residues between WT and CUA share strictly identical coordinates ($\text{RMSD} = 0.000\ \text{Å}$), eliminating structural confounding artifacts in downstream comparative MD.
- **Documentation:** Full methodology recorded in [docs/12_structure_modeling_method.md](file:///d:/trna-asl-md-stability/docs/12_structure_modeling_method.md).

---

### SDR-12: Pre-MD Structural & Force Field Validation Standards
- **Status:** APPROVED & IMPLEMENTED `[ENGINEERING CHECK]`
- **Context:** Validating starting coordinates before simulation requires distinguishing physical truths from operational engineering thresholds.
- **Decision:** Establish a 4-tier validation suite in `scripts/02_validate_asl_structure.py`:
  1. *Tier 1 (Biological Integrity `[HARD REQUIREMENT]`):* 17-nt length, Sprinzl numbering 27–43, canonical RNA residues, sequence match, invariant U33, purine A37.
  2. *Tier 2 (Geometric Sanity `[ENGINEERING CHECK]`):* Bond lengths within standard covalent bounds, glycosidic $\text{anti}$ torsions, $C3'$-endo sugar puckers, absence of steric clashes ($> 1.80\ \text{Å}$ global, $> 3.2\ \text{Å}$ local mutation site clearance).
  3. *Tier 3 (Comparative Control `[SCIENTIFIC DECISION]`):* Mutations strictly confined to 34 and 36; all 15 non-mutated residues have $\text{RMSD} = 0.000000\ \text{Å}$.
  4. *Tier 4 (Force Field Parameterization `[ENGINEERING CHECK]`):* OpenMM Modeller hydrogen resolution (544 WT / 543 CUA atoms) and AMBER OL3 parameterization (`createSystem()`) verifying exact net charge of $-16.0000\ e$.
- **Explicit Boundary:** Pre-MD validation proves initial coordinate validity only; it does NOT prove dynamic stability or predict suppression efficiency.

---

### SDR-13: Two-Stage Solvated Energy Minimization Protocol
- **Status:** APPROVED & IMPLEMENTED (PHASE 4B) `[IMPLEMENTATION CHOICE]`
- **Context:** Solvating rigid starting structures introduces high-energy steric overlaps and unrelaxed electrostatic contacts between solvent molecules, ions, and the RNA solute. Sudden unrestrained minimization of an isolated RNA hairpin in explicit solvent risks premature deformation of the native A-form stem or crystallographic U-turn loop before the solvent sheath relaxes.
- **Decision:** Implement a two-stage minimization protocol in `scripts/05_minimize.py`:
  1. *Stage 1 (Restrained Solute Minimization):* Apply harmonic position restraints ($k = 5.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2} \approx 2092.00\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$) to all RNA heavy atoms while allowing TIP3P water molecules and $\text{K}^+/\text{Cl}^-$ ions to minimize (max 2,500 iterations, tolerance $10.0\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$).
  2. *Stage 2 (Full System Minimization):* Release position restraints ($k = 0.0$) and minimize the entire solvated system (max 5,000 iterations, tolerance $10.0\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$).
- **Rationale & Empirical Validation:**
  - Stage 1 accounts for $>99\%$ of the total potential energy drop ($\approx -167,500\ \text{kJ/mol}$ out of $\approx -169,000\ \text{kJ/mol}$), relaxing solvent packing without distorting the solute.
  - Stage 2 relieves remaining local steric strain within the mutated nucleobases ($C34$ and $A36$) in the context of relaxed solvent.
  - RNA heavy-atom RMSD relative to starting coordinates is preserved under $1.0\ \text{Å}$ ($0.6346\ \text{Å}$ WT, $0.8921\ \text{Å}$ CUA), with zero steric clashes ($> 1.22\ \text{Å}$).
- **Explicit Boundary:** Energy minimization removes unfavorable local geometries but does not establish conformational stability, folding, or biological function. Full dynamic evaluation is deferred to production MD.

---

### SDR-14: Staged NVT/NPT Restrained Equilibration Protocol
- **Status:** APPROVED & IMPLEMENTED (PHASE 4C) `[IMPLEMENTATION CHOICE]`
- **Context:** Following energy minimization, the system is at 0 K with frozen solvent packing. Sudden unconstrained dynamics at 310 K and 1 bar would introduce thermal shock, solvent cavitation, and potentially unphysical distortion of the RNA hairpin before solvent density and box volume relax to equilibrium values.
- **Decision:** Implement a multi-stage restrained equilibration pipeline in `scripts/06_equilibrate.py`:
  1. *Part I: NVT Restrained Heating (10 K -> 310 K, 200 ps):* Apply harmonic restraints ($k = 5.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$) to all 357/356 RNA heavy atoms. Initialize velocities from Maxwell-Boltzmann distribution at 10.0 K (deterministic seeds: WT = 4101, CUA = 4201). Heat linearly to 310.0 K over 150 ps (75,000 steps at $\Delta t = 2.0\ \text{fs}$) and plateau at 310.0 K for 50 ps using `LangevinMiddleIntegrator` ($\gamma = 1.0\ \text{ps}^{-1}$).
  2. *Part II: NPT Density Relaxation & Progressive Restraint Release (310 K, 1.0 bar, 1000 ps):* Couple to OpenMM `MonteCarloBarostat` (attempted every 25 steps) across 4 sequential stages (250 ps each):
     - Stage 1: $k = 2.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ (rapid box condensation from initial $166.4\ \text{nm}^3$ to $\sim 159.4\ \text{nm}^3$).
     - Stage 2: $k = 1.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ (solvent density relaxation).
     - Stage 3: $k = 0.2\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ (mild guiding restraints allowing backbone relaxation).
     - Stage 4: $k = 0.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ (**genuinely unrestrained pre-production relaxation**).
- **Rationale & Empirical Validation:**
  - Total equilibration time is 1.2 ns (600,000 integration steps per construct).
  - Bulk solvent density converges cleanly to $1.0182\text{--}1.0193\ \text{g/cm}^3$ (physiologically correct density of 150 mM KCl water at 310 K).
  - RNA heavy-atom Kabsch RMSD relative to minimized starting structures is preserved: WT = 1.66 Å (stem: 0.88 Å), CUA = 1.92 Å (stem: 0.93 Å), demonstrating that the A-form helical stem and loop architecture remain intact with zero steric clashes.
- **Explicit Boundary:** Equilibration relaxes solvent density and heats the system to target thermodynamic conditions, but does NOT establish conformational stability, native folding, or biological function. Full dynamic stability and comparative statistical testing are deferred to Phase 5 production MD.
