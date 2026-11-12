# 04 - Scientific Decision Log

This log formally records all scientific and technical decisions, their rationales, alternatives considered, and governance status. In accordance with project instructions, decisions mandated by the take-home are distinguished from implementation choices and open decisions.

---

## Record Index

| ID | Title | Status | Classification |
| :--- | :--- | :--- | :--- |
| **SDR-01** | ASL Structural Boundary (17-nt, 5-bp Stem, 7-nt Loop) | **APPROVED** | `[PROMPT MANDATE]` |
| **SDR-02** | Force Field Parameterization: AMBER OL3 (`ff99bsc0chiOL3`) | **APPROVED** | `[PROMPT MANDATE]` |
| **SDR-03** | Thermodynamic Ensemble (NPT at 310 K, 1 bar) | **APPROVED** | `[PROMPT MANDATE]` |
| **SDR-04** | Sampling Budget (3 replicates × 20 ns per construct = 120 ns) | **APPROVED** | `[PROMPT MANDATE]` |
| **SDR-05** | Ionic Environment (~0.15 M KCl in Explicit Water) | **APPROVED** | `[PROMPT MANDATE]` |
| **SDR-06** | Solvent Model & Matched Ion Parameter Selection | **BASELINE ESTABLISHED** | `[IMPLEMENTATION CHOICE]` |
| **SDR-07** | Target Human Cytoplasmic tRNA Selection Criteria | **DEFERRED TO PHASE 2** | `[PENDING VERIFICATION]` |
| **SDR-08** | Selection & Justification of Metric 6 | **DEFERRED TO PHASE 2** | `[PENDING VERIFICATION]` |
| **SDR-09** | ASL Terminal Chemistry (5'-OH / 3'-OH) | **PROPOSED BASELINE** | `[SCIENTIFIC ASSUMPTION]` |
| **SDR-10** | Replicate-Aware Statistical Inference & Sample Size Constraints | **METHODOLOGY DEFINED** | `[SCIENTIFIC ASSUMPTION]` |

---

## Detailed Scientific Decision Records

### SDR-01: ASL Structural Boundary Definition
- **Status:** APPROVED `[PROMPT MANDATE]`
- **Context:** An anticodon stem-loop (ASL) within a canonical tRNA spans residues 27–43 in Sprinzl numbering.
- **Decision:** Model an isolated 17-nucleotide hairpin comprising a 5-bp Watson-Crick stem (base pairs 27:43, 28:42, 29:41, 30:40, 31:39) and a 7-nucleotide anticodon loop (residues 32–38).
- **Rationale:** Minimal autonomous structural unit representing the anticodon domain mandated by the take-home requirements.

---

### SDR-02: Force Field Parameterization (AMBER OL3)
- **Status:** APPROVED `[PROMPT MANDATE]`
- **Context:** Empirical force fields require accurate glycosidic torsions to avoid unphysical RNA laddering or loop collapse.
- **Decision:** Utilize the AMBER OL3 parameter set (`ff99bsc0chiOL3`).
- **Rationale:** OL3 specifically corrects RNA $\chi$ glycosidic torsion biases and represents the established standard for RNA hairpins.

---

### SDR-03: Thermodynamic Ensemble (NPT at 310 K, 1 bar)
- **Status:** APPROVED `[PROMPT MANDATE]`
- **Context:** Human cytoplasmic tRNAs function in vivo at physiological body temperature.
- **Decision:** Simulate in the isothermal-isobaric ensemble ($NPT$) at $T = 310.0\ \text{K}$ and $P = 1.0\ \text{bar}$.
- **Rationale:** Matches physiological human body temperature ($37^\circ\text{C}$) and standard pressure.

---

### SDR-04: Sampling Budget & Replicate Control
- **Status:** APPROVED `[PROMPT MANDATE]`
- **Context:** RNA hairpin trajectories can become trapped in local kinetic minima.
- **Decision:** Execute 3 statistically independent replicates per construct (WT and CUA suppressor), each running for a minimum of 20 ns production, totaling 6 simulations and 120 ns aggregate sampling. Assign initial velocities using distinct pseudo-random seeds.
- **Rationale:** Multiple independent trajectories provide an assessment of replicate variability that cannot be captured by a single trajectory.

---

### SDR-05: Ionic Environment (~0.15 M KCl)
- **Status:** APPROVED `[PROMPT MANDATE]`
- **Context:** RNA is a polyanion requiring monovalent counterions to screen phosphate repulsions.
- **Decision:** Solvate in explicit water, neutralize the 16 negative charges with $\text{K}^+$, and add $\text{K}^+$ and $\text{Cl}^-$ co-ions to reach an ionic strength of approximately $0.15\ \text{M}$.
- **Rationale:** Matches physiological human cytoplasmic potassium concentrations ($\sim 140\text{--}150\ \text{mM}\ \text{K}^+$).

---

### SDR-06: Solvent Model & Matched Ion Parameter Selection
- **Status:** BASELINE ESTABLISHED `[IMPLEMENTATION CHOICE]`
- **Context:** The take-home prompt explicitly permits TIP3P or an equivalent model. The water model must be compatible with AMBER OL3 and matched with appropriate ion parameters.
- **Decision:** 
  - Establish **TIP3P** as the primary baseline water model for this project because AMBER OL3 (`ff99bsc0chiOL3`) was parameterized and validated natively in TIP3P.
  - Pair TIP3P with matched monovalent ion parameters (e.g., Joung-Cheatham `ionsjc_tip3p` or standard AMBER ion libraries).
  - Alternative water models (such as OPC or TIP4P-EW) remain optional only if accompanied by explicit, documented scientific justification and strictly matched ion parameters (e.g., `ionsjc_opc`).
- **Rationale:** Enforces internal consistency across force field, solvent, and ion parameters, avoiding arbitrary or unverified solvent preferences.

---

### SDR-07: Target Human Cytoplasmic tRNA Selection Criteria
- **Status:** DEFERRED TO PHASE 2 `[PENDING VERIFICATION]`
- **Context:** The prompt requires comparing a wild-type human cytoplasmic tRNA ASL against its CUA suppressor mutant, under the strict constraint that the native anticodon must convert to `CUA` by **exactly two nucleotide substitutions** ($\text{Hamming}(\text{Anticodon}_{\text{WT}}, \text{CUA}) = 2$).
- **Audit Retraction:** 
  - Previously suggested candidates (tRNA-Tyr, tRNA-Trp, tRNA-Leu) in the initial scaffold draft are **retracted** because:
    1. tRNA-Tyr (`GUA` $\rightarrow$ `CUA`) requires only 1 mutation ($G34C$).
    2. tRNA-Trp (`CCA` $\rightarrow$ `CUA`) requires only 1 mutation ($C35U$).
    3. tRNA-Leu (`CAA` $\rightarrow$ `CUA`) requires only 1 mutation ($A35U$).
    4. Previous informal PDB citations were misattributed non-human structures and have been completely removed.
- **Phase 2 Selection Protocol:**
  - No tRNA is selected in Phase 1.
  - In Phase 2, candidate genes will be screened against verified public records in GtRNAdb for human cytoplasmic tRNAs.
  - Screening filters:
    1. Canonical 5-bp anticodon stem (residues 27–31 paired with 39–43).
    2. Canonical 7-nt loop (residues 32–38) with invariant U33.
    3. Wild-type anticodon triplet exhibiting $\text{Hamming}(\text{Anticodon}_{\text{WT}}, \text{CUA}) = 2$.
    4. Availability of experimental structural templates or high-confidence homology models.

---

### SDR-08: Selection & Justification of Metric 6
- **Status:** DEFERRED TO PHASE 2 `[PENDING VERIFICATION]`
- **Context:** The take-home prompt requires one additional justified structural metric alongside RMSD, RMSF, H-bonds, U-turn diagnostic, and 3'-stacking.
- **Decision:** No metric is prematurely locked in Phase 1. Formal selection and scientific justification will occur in Phase 2 from the following candidates:
  - *Candidate A:* **Anticodon Loop Radius of Gyration ($R_{g,\text{loop}}$)** — Quantifies global loop expansion or compaction.
  - *Candidate B:* **Anticodon Solvent-Accessible Surface Area (SASA)** — Quantifies nucleobase solvent exposure vs. lumen flipping.
  - *Candidate C:* **Stem-Loop Interhelical Angle / Axis Tilt** — Quantifies coaxial bending relative to the stem helix.

---

### SDR-09: ASL Terminal Chemistry
- **Status:** PROPOSED BASELINE `[SCIENTIFIC ASSUMPTION]`
- **Context:** Simulating an excised 17-nt hairpin requires defining terminal groups at residues 27 and 43.
- **Decision:** Maintain neutral **5'-OH** and **3'-OH** termini as the proposed baseline.
- **Rationale:** Prevents placing a repulsive $-2e$ charge on residue 27 directly opposite residue 43, minimizing artificial terminal base-pair fraying. Final confirmation will occur during Phase 2 coordinate assembly.

---

### SDR-10: Replicate-Aware Statistical Inference & Sample Size Constraints
- **Status:** METHODOLOGY DEFINED `[SCIENTIFIC ASSUMPTION]`
- **Context:** Trajectory frames sampled at 10 ps intervals exhibit high temporal autocorrelation. Furthermore, with $N=3$ independent replicates per construct, statistical power for traditional null-hypothesis testing is severely constrained.
- **Mathematical Correction from Audit:**
  - A two-sided Mann-Whitney $U$ test with $N_1 = 3, N_2 = 3$ has $\binom{6}{3} = 20$ possible rank arrangements. The minimum attainable two-sided $p$-value is $2/20 = 0.10$.
  - Therefore, claiming that classical significance testing can achieve $p < 0.05$ with $N=3$ is mathematically impossible.
- **Decision:**
  1. **Primary Reporting:** Center reporting on replicate-level effect sizes (Hedges' $g$ / Cohen's $d$, difference in replicate means/medians) with 95% bootstrap confidence intervals ($10,000$ iterations).
  2. **Explicit Disclosure:** Transparently state sample-size and power limitations in the final report.
  3. **Pseudoreplication Guard:** Trajectory frames will not be treated as independent biological replicates. Autocorrelation times ($\tau_{\text{int}}$) and block averages will be computed to assess intra-trajectory stability, while inter-replicate comparisons remain the primary statistical unit.
