# 04 - Scientific Decision Log

This log formally records all scientific and technical decisions, their rationales, alternatives considered, and governance status. In accordance with project instructions, no scientific decision is made silently.

---

## Record Index

| ID | Title | Status | Scope |
| :--- | :--- | :--- | :--- |
| **SDR-01** | ASL Structural Boundary (17-nt, 5-bp Stem, 7-nt Loop) | **APPROVED** (Prompt) | System Architecture |
| **SDR-02** | Force Field Parameterization: AMBER OL3 (`ff99bsc0chiOL3`) | **APPROVED** (Prompt) | Molecular Mechanics |
| **SDR-03** | Thermodynamic Ensemble (NPT at 310 K, 1 bar) | **APPROVED** (Prompt) | Simulation Protocol |
| **SDR-04** | Sampling Budget (3 replicates × 20 ns per construct = 120 ns) | **APPROVED** (Prompt) | Simulation Protocol |
| **SDR-05** | Ionic Environment (~0.15 M KCl in Explicit Water) | **APPROVED** (Prompt) | Solvation & Physics |
| **SDR-06** | Water Model Selection (OPC vs. TIP3P vs. TIP4P-EW) | **PENDING USER APPROVAL** | Solvation & Physics |
| **SDR-07** | Selection of Candidate Human Cytoplasmic tRNA | **PENDING USER APPROVAL** | Biological Construct |
| **SDR-08** | Selection & Justification of Metric 6 | **PENDING USER APPROVAL** | Trajectory Analysis |
| **SDR-09** | ASL Terminal Chemistry (5'-OH / 3'-OH vs. 5'-monophosphate) | **PENDING USER APPROVAL** | Structure Preparation |
| **SDR-10** | Statistical Autocorrelation & Replicate Pooling Methodology | **PENDING USER APPROVAL** | Statistical Inference |

---

## Detailed Scientific Decision Records

### SDR-01: ASL Structural Boundary Definition
- **Status:** APPROVED (Mandated by project specification)
- **Context:** An anticodon stem-loop (ASL) within a full tRNA spans residues 27–43 in canonical Sprinzl numbering.
- **Decision:** Model an isolated 17-nucleotide hairpin comprising a 5-bp Watson-Crick stem (base pairs 27:43, 28:42, 29:41, 30:40, 31:39) and a 7-nucleotide anticodon loop (residues 32–38).
- **Rationale:** The isolated 17-nt ASL is the established minimal autonomous folding domain for studying anticodon dynamics and U-turn stability in biophysical literature.
- **Alternatives Considered:** 21-nt ASL (including base pairs 26:44). Rejected because prompt mandates a 5-bp stem and 17-nt length.

---

### SDR-02: Force Field Parameterization (AMBER OL3)
- **Status:** APPROVED (Mandated by project specification)
- **Context:** Classical force fields (e.g., standard AMBER ff99 or CHARMM27) suffer from known artifacts in RNA simulations, such as non-canonical $\alpha/\zeta$ and $\chi$ torsion flips leading to unphysical laddering or loop collapse.
- **Decision:** Utilize the AMBER OL3 parameter set (`ff99bsc0chiOL3`).
- **Rationale:** OL3 provides re-parameterized glycosidic $\chi$ torsions that accurately maintain RNA loop geometries and A-form helices, and is recognized as a benchmark standard in nucleic acid MD.
- **Alternatives Considered:** CHARMM36, DES-Amber. AMBER OL3 is explicitly mandated by the project requirements.

---

### SDR-03: Thermodynamic Ensemble (NPT at 310 K, 1 bar)
- **Status:** APPROVED (Mandated by project specification)
- **Context:** Human cytoplasmic tRNAs function in vivo at physiological body temperature.
- **Decision:** Simulate in the isothermal-isobaric ensemble ($NPT$) at $T = 310.0\ \text{K}$ and $P = 1.0\ \text{bar}$.
- **Rationale:** Matches physiological human body temperature ($37^\circ\text{C}$) and physiological atmospheric pressure, ensuring biological relevance of conformational sampling.

---

### SDR-04: Sampling Budget & Replicate Control
- **Status:** APPROVED (Mandated by project specification)
- **Context:** RNA hairpin trajectories can become trapped in local kinetic minima.
- **Decision:** Execute 3 statistically independent replicates per construct (WT and CUA suppressor), each running for a minimum of 20 ns production, totaling 6 simulations and 120 ns aggregate sampling. Assign initial velocities using distinct pseudo-random seeds.
- **Rationale:** Multiple independent trajectories provide an assessment of conformational divergence and replicate reproducibility that cannot be captured by a single long trajectory.

---

### SDR-05: Ionic Environment (~0.15 M KCl)
- **Status:** APPROVED (Mandated by project specification)
- **Context:** RNA is a highly charged polyanion requiring monovalent counterions to screen phosphate repulsions.
- **Decision:** Solvate with explicit water and neutralize the 16 negative phosphate charges with $\text{K}^+$, then add $\text{K}^+$ and $\text{Cl}^-$ co-ions to reach an ionic strength of approximately $0.15\ \text{M}$.
- **Rationale:** Matches intracellular human cytoplasmic potassium concentrations ($\sim 140\text{--}150\ \text{mM}\ \text{K}^+$), reflecting physiological conditions superior to unphysiological $\text{NaCl}$.

---

### SDR-06: Water Model Selection
- **Status:** PENDING USER APPROVAL
- **Context:** Water models strongly modulate nucleic acid conformational equilibria. Standard TIP3P is known to overestimate RNA base-stacking and allow artificial loop compaction over longer timescales. Four-point (TIP4P-EW, TIP4P-D) or optimal four-point (OPC) models improve RNA conformational sampling.
- **Proposed Options:**
  - *Option A (Recommended):* **OPC** (Optimal Point Charge). Shown in literature (e.g., Shabane et al., J. Chem. Theory Comput. 2019) to yield superior RNA hairpin and loop geometries when paired with AMBER OL3.
  - *Option B:* **TIP4P-EW**. Well-balanced electrostatic screening, widely validated with AMBER nucleic acid parameters.
  - *Option C:* **TIP3P**. Traditional, highly compatible, but higher tendency toward artificial loop compaction.
- **Action Required:** User approval on preferred solvent model prior to Phase 3 system setup.

---

### SDR-07: Selection of Candidate Wild-Type Human Cytoplasmic tRNA
- **Status:** PENDING USER APPROVAL
- **Context:** The prompt requires comparing a wild-type human cytoplasmic tRNA ASL against its CUA suppressor mutant. No specific tRNA was pre-selected.
- **Proposed Candidates (Public Data Only):**
  - *Candidate 1 (Recommended):* **Human cytoplasmic tRNA-Tyr (GTA / UUA)**
    - High historical relevance: tRNA-Tyr is the archetypal amber suppressor studied extensively in genetic code expansion (engineered amber suppressor tRNAs).
    - Clear 5-bp stem and canonical U33-containing loop.
    - Abundant high-resolution crystallographic reference structures available in PDB (e.g., PDB ID: 4TRA, 1EHZ, 6U8D homologs).
  - *Candidate 2:* **Human cytoplasmic tRNA-Trp (CCA)**
    - Single mutation at position 34 ($C \rightarrow U$) converts the native `CCA` anticodon to `UCA` (opal) or dual mutation ($C34U, C35A$) yields `CUA` (amber).
  - *Candidate 3:* **Human cytoplasmic tRNA-Leu (CAA)**
    - Well-characterized human cytoplasmic tRNA; suppression efficiency studied in synthetic biology.
- **Action Required:** User decision on candidate tRNA selection prior to Phase 2 modeling.

---

### SDR-08: Selection & Justification of Metric 6
- **Status:** PENDING USER APPROVAL
- **Context:** The take-home prompt requires one additional justified structural metric alongside RMSD, RMSF, H-bonds, U-turn geometry, and 3'-stacking.
- **Proposed Candidates:**
  - *Candidate A (Recommended):* **Anticodon Loop Radius of Gyration ($R_{g,\text{loop}}$)**
    - *Scientific Justification:* Directly measures overall loop compactness vs. expansion. Destabilization of the canonical U-turn motif or loss of 3'-stacking invariably leads to loop swelling (increased $R_g$) or unphysical collapse (compacted $R_g$).
  - *Candidate B:* **Anticodon Solvent-Accessible Surface Area (SASA)**
    - *Scientific Justification:* Measures exposure of nucleobases 34, 35, and 36 to solvent, reflecting whether bases remain projected outward (as required for codon reading) or flip inward into the loop lumen.
  - *Candidate C:* **Stem-Loop Interhelical Angle / Axis Tilt**
    - *Scientific Justification:* Evaluates whether the anticodon loop maintains coaxial alignment with the stem helix or tilts abnormally upon mutation.
- **Action Required:** User approval of Metric 6 selection prior to Phase 5 analysis.

---

### SDR-09: ASL Terminal Chemistry
- **Status:** PENDING USER APPROVAL
- **Context:** When simulating an isolated 17-nt RNA fragment excised from a full tRNA, the terminal residues (27 and 43) must be capped or terminated.
- **Proposed Options:**
  - *Option A (Recommended):* Neutral **5'-OH** and **3'-OH** termini.
    - *Scientific Justification:* Avoids placing a full negative charge (5'-monophosphate) at residue 27 immediately adjacent to residue 43, preventing artificial electrostatic fraying of the closing base pair 27:43.
  - *Option B:* Standard **5'-monophosphate** and **3'-OH**.
    - *Scientific Justification:* Reflects RNA cleavage products, but adds $-2e$ charge to terminal nucleotide.
- **Action Required:** User confirmation of terminal chemistry.

---

### SDR-10: Statistical Autocorrelation & Replicate Pooling
- **Status:** PENDING USER APPROVAL
- **Context:** Molecular dynamics frames sampled every 10 ps exhibit significant temporal autocorrelation. Treating frames as independent artificially inflates degrees of freedom ($p$-hacking artifact).
- **Proposed Methodology:**
  1. Compute statistical inefficiency / autocorrelation time ($\tau_{\text{int}}$) for each metric using the automated windowing algorithm (Chodera et al.).
  2. Subsample trajectories at intervals $\Delta t \ge 2\tau_{\text{int}}$ or apply block averaging to calculate effective independent samples ($N_{\text{eff}}$).
  3. Treat each of the 3 replicates as the primary statistical unit (replicate-level means), computing standard error across replicates and testing hypothesis differences using non-parametric permutation tests or Mann-Whitney $U$ tests.
- **Action Required:** User confirmation of statistical inference methodology.
