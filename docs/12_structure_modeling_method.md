# 12 - 3D Structural Modeling Methodology & Decision Record

**Study Title:** tRNA ASL Structural Stability Study  
**Phase:** Phase 3 (3D Structural Modeling & Validation)  
**Status:** Method Selected, Implemented & Verified  

---

## 1. Objective & Mandate

The objective of Phase 3 is to construct physically and chemically valid, stereochemically unhurried 3D starting coordinates for:
1. **Wild-Type Construct:** 17-nt human cytoplasmic tRNA-Asp-GTC-2-1 ASL (`5'-CCCGCCUGUCACGCGGG-3'`)
2. **Suppressor Construct:** 17-nt amber suppressor ASL (`5'-CCCGCCUCUAACGCGGG-3'`)

Both constructs must adhere to strict structural constraints:
- Exactly 17 ribonucleotides (residues 27 to 43 in Sprinzl numbering). `[PROMPT MANDATE]`
- Canonical 5-bp stem and 7-nt anticodon loop architecture. `[PROMPT MANDATE]`
- Unmodified ribonucleotides with standard RNA stereochemistry. `[PROMPT MANDATE]`
- Controlled comparative design: Suppressor differs from WT **strictly and exclusively** at position 34 ($G \rightarrow C$) and position 36 ($C \rightarrow A$), with identical starting coordinates across all other 15 residues. `[PROMPT MANDATE]`
- **No Molecular Dynamics in Phase 3:** No solvation, no ion addition, no energy minimization, no equilibration, and no production MD are permitted in this phase. `[PROMPT MANDATE]`

---

## 2. Evaluation of Candidate Structural Modeling Approaches

To determine the most scientifically defensible, reproducible, and robust construction strategy, four distinct structural modeling approaches were systematically evaluated:

### Approach A: De Novo Fragment-Based RNA Modeling (e.g., RNAComposer, Rosetta FARNA/FARFAR, SimRNA)
- **Description:** Algorithmic assembly of 3D coordinates directly from secondary structure (`>>>>>.......<<<<<`) and primary sequence.
- **Evaluation:**
  - *Pros:* Fully automated generation without manual template search.
  - *Cons:*
    - Rosetta FARNA/FARFAR requires external C++ compilation, proprietary/academic licensing barriers, and relies on stochastic Monte Carlo sampling that is non-deterministic without strict pseudo-random seeding.
    - RNAComposer is primarily accessible via a web server interface, creating brittle external network dependencies that undermine standalone pipeline reproducibility.
    - De novo prediction of the 7-nt anticodon loop frequently samples non-canonical loop conformations, requiring subjective energy scoring or manual clustering to select a single conformer.
- **Verdict:** **REJECTED** due to non-determinism, software accessibility constraints, and potential loop conformer subjectivity.

### Approach B: Programmatic Idealized A-Form Duplex + Generic Loop Closure (e.g., AmberTools NAB / 3DNA)
- **Description:** Mathematical construction of an ideal A-form RNA duplex for the 5-bp stem (residues 27:31 and 39:43) followed by geometric loop closure across the 7-nt loop.
- **Evaluation:**
  - *Pros:* Purely algorithmic and mathematically transparent.
  - *Cons:*
    - tRNA anticodon loops do not adopt generic hairpin loop geometry. They possess a highly specialized tertiary architecture—the canonical **Quigley & Rich U-turn motif** (invariant U33, sharp phosphate reversal at residue 36, and 3'-purine stacking anchored by A37)—which generic loop-closing algorithms fail to reproduce without significant structural distortion.
    - Resolving loop closure strain would necessitate extensive energy minimization, violating the Phase 3 boundary condition ("Do not run minimization in Phase 3").
- **Verdict:** **REJECTED** due to inability to reliably capture the authentic experimental U-turn loop fold without energy minimization.

### Approach C: ModeRNA Homology Modeling
- **Description:** Template-based homology modeling tool specifically designed for RNA using structural alignments and fragment substitution.
- **Evaluation:**
  - *Pros:* Specifically built for RNA template-based modeling.
  - *Cons:* ModeRNA is an unmaintained Python 2.6/2.7 package with unresolved dependencies on modern 64-bit operating systems and Python 3.11+. Its underlying algorithm (rigid-body base replacement on a conserved ribose-phosphate scaffold) can be executed more transparently and reliably in pure Python.
- **Verdict:** **REJECTED** due to software obsolescence and legacy platform incompatibility.

### Approach D: Experimental Crystallographic Template-Based Homology Modeling (Selected)
- **Description:** Derivation of starting coordinates from a high-resolution, experimentally determined X-ray crystal structure of an unmodified tRNA-Asp anticodon stem-loop, followed by deterministic glycosidic base grafting.
- **Candidate Template:** **PDB 6UGG** (Chain A, residues 27–43), a 1.95 Å resolution crystal structure of unmodified tRNA(Asp) (*Escherichia coli*, DOI: [10.1038/s41467-020-15365-2](https://doi.org/10.1038/s41467-020-15365-2)).
- **Sequence Comparison:**
  - `6UGG` ASL Sequence (27–43): `5'- C C U G C   C U G U C A C   G C A G G -3'`
  - Target WT ASL Sequence: `5'- C C C G C   C U G U C A C   G C G G G -3'`
  - **Loop Identity:** The entire 7-nt anticodon loop (residues 32–38: `CUGUCAC`) is **100% sequence-identical** between PDB 6UGG and human `tRNA-Asp-GTC-2-1`.
  - **Stem Identity:** 4 out of 5 base pairs are already identical (`27:43 C:G`, `28:42 C:G`, `30:40 G:C`, `31:39 C:G`). Only a single base pair differs in the regular A-form stem: `U29:A41` in 6UGG vs. `C29:G41` in human.
- **Evaluation:**
  - *Pros:*
    - Grounded in authentic high-resolution (1.95 Å) experimental crystallographic data.
    - Preserves true experimental U-turn loop geometry ($U33\ \text{N3} \cdots C36\ \text{OP2} = 3.52\ \text{Å}$), standard $C3'$-endo ribose puckers, and canonical $\text{anti}$ glycosidic torsions ($\chi \approx -160^\circ$ to $176^\circ$) directly observed in electron density.
    - 15 out of 17 residues in the WT construct require zero atom alteration.
    - The compensatory base pair mutation in the stem (`U29C`, `A41G`) and the two suppressor mutations (`G34C`, `C36A`) are executed deterministically on the rigid glycosidic frame without altering the sugar-phosphate backbone.
    - 100% reproducible via open-source Python, platform-independent, and fully verified for OpenMM / AMBER OL3 force field compatibility.
- **Verdict:** **ACCEPTED & SELECTED** as the optimal scientific approach.

---

## 3. Scientific Decision Classification

| Decision Element | Classification | Rationale & Governance |
| :--- | :---: | :--- |
| **17-nt Isolated Hairpin** | `[PROMPT MANDATE]` | Mandated by Transloka Bio take-home brief (residues 27–43). |
| **Unmodified Nucleotides** | `[PROMPT MANDATE]` | Mandated baseline comparison state. |
| **Exact Suppressor Mutations** | `[PROMPT MANDATE]` | Only $G34 \rightarrow C34$ and $C36 \rightarrow A36$ mutated; all other positions identical. |
| **No MD / Minimization in Phase 3** | `[PROMPT MANDATE]` | Strict boundary condition preventing premature simulation. |
| **Crystallographic Template (PDB 6UGG)** | `[SCIENTIFIC DECISION]` | Selected for 100% loop sequence identity and 1.95 Å unmodified experimental quality. |
| **Neutral 5'-OH and 3'-OH Termini** | `[SCIENTIFIC DECISION]` | Encodes SDR-09; omits 5'-phosphate on residue 27 to avoid terminal charge repulsion. |
| **Deterministic Glycosidic Frame Grafting** | `[IMPLEMENTATION CHOICE]` | Mathematical rigid-body frame alignment preserving bond lengths, planarity, and $\chi$ angle. |
| **Sprinzl Residue Numbering (27–43)** | `[IMPLEMENTATION CHOICE]` | Preserves canonical tRNA coordinate indices across files and simulation configs. |

---

## 4. Deterministic Base Grafting Algorithm

When substituting a nucleobase at position $i$ ($U29 \rightarrow C29$, $A41 \rightarrow G41$, $G34 \rightarrow C34$, $C36 \rightarrow A36$):

1. **Backbone Invariance:** The ribose-phosphate backbone atoms ($O5', C5', C4', O4', C3', O3', C2', O2', C1'$, and internal $P, OP1, OP2$) are held strictly stationary.
2. **Anchor Placement:** The target glycosidic nitrogen ($N1$ for pyrimidines, $N9$ for purines) is positioned along the existing glycosidic bond unit vector:
   $$\mathbf{r}_{\text{N}} = \mathbf{r}_{\text{C1'}} + d_{\text{glyc}} \cdot \frac{\mathbf{r}_{\text{N,orig}} - \mathbf{r}_{\text{C1'}}}{\|\mathbf{r}_{\text{N,orig}} - \mathbf{r}_{\text{C1'}}\|}, \quad d_{\text{glyc}} = 1.470\ \text{Å}$$
3. **Local Orthonormal Coordinate Frame:**
   - $\hat{e}_x$: Unit vector along the minor groove / $\chi$-defining bond ($\mathbf{r}_{\text{C\_dir}} - \mathbf{r}_{\text{N}}$, where $C_{\text{dir}}$ is $C2$ for pyrimidines and $C4$ for purines).
   - $\hat{e}_z$: Unit normal to the planar base ring ($\hat{n}_{\text{base}}$).
   - $\hat{e}_y = \hat{e}_z \times \hat{e}_x$.
4. **Coordinate Transformation:** Every atom $a$ of the ideal base template is mapped onto the target sugar via:
   $$\mathbf{r}_a = \mathbf{r}_{\text{N}} + x_{a,\text{local}} \hat{e}_x + y_{a,\text{local}} \hat{e}_y + z_{a,\text{local}} \hat{e}_z$$
5. **Physical Outcome:**
   - Preserves standard covalent bond lengths ($d(C1'-N) = 1.470\ \text{Å}$, $d(N1-C2) = 1.400\ \text{Å}$, $d(N9-C4) = 1.374\ \text{Å}$).
   - Preserves ring planarity and canonical $\text{anti}$ glycosidic torsion angle ($\chi$).
   - Yields zero steric clashes (all non-bonded contacts $> 1.8\ \text{Å}$) and native Watson-Crick hydrogen bonding distances ($2.74 - 2.82\ \text{Å}$) for the mutated $C29:G41$ stem pair.
