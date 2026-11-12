# 01 - Project Specification

## 1. Project Title & Take-Home Research Question

**Project Title:** tRNA ASL Structural Stability Study  
**Take-Home Research Question:** `[PROMPT MANDATE]`
> *"When the anticodon of a human cytoplasmic tRNA is mutated to a suppressor anticodon, does the anticodon stem-loop (ASL) retain canonical loop architecture or become destabilized?"*

This study investigates whether introducing an amber suppressor mutation (`CUA` anticodon, targeting the `UAG` stop codon) into the anticodon stem-loop (ASL) of a human cytoplasmic tRNA alters its structural dynamics, stem integrity, or canonical loop architecture relative to the wild-type (WT) sequence.

---

## 2. Comparative Experimental Design

The study is formulated as a controlled, comparative computational experiment between two molecular constructs:

1. **Construct 1 (Wild-Type ASL):**
   - The wild-type 17-nucleotide ASL sequence derived from verified human cytoplasmic **tRNA-Asp-GTC-2-1** in GtRNAdb (`data/metadata/selected_trna.yaml`, `docs/11_trna_selection.md`). `[VERIFIED]`
   - **Anticodon Conversion Requirement:** `[PROMPT MANDATE]` The wild-type anticodon triplet (positions 34, 35, 36) must differ from `CUA` by **exactly two nucleotide substitutions** ($\text{Hamming}(\text{Anticodon}_{\text{WT}}, \text{CUA}) = 2$).
2. **Construct 2 (Suppressor ASL):**
   - An identical 17-nucleotide sequence except for the replacement of the native anticodon triplet with the `CUA` suppressor anticodon at positions 34–36. `[PROMPT MANDATE]`
   - The stem base pairs (residues 27–31 and 39–43) and loop flanking positions (residues 32, 33, 37, 38) are preserved identically to isolate the specific structural impact of the anticodon mutation. `[SCIENTIFIC ASSUMPTION]`

---

## 3. Structural Model Specifications

| Parameter | Specification | Classification | Scientific Rationale |
| :--- | :--- | :--- | :--- |
| **System Size** | 17 nucleotides | `[PROMPT MANDATE]` | Minimal autonomous structural unit representing the tRNA anticodon domain. |
| **Stem Architecture** | 5 base pairs | `[PROMPT MANDATE]` | Canonical secondary structure of human cytoplasmic tRNA anticodon stems (residues 27–31 paired with 39–43 in Sprinzl numbering). |
| **Loop Architecture** | 7 nucleotides | `[PROMPT MANDATE]` | Canonical heptanucleotide loop (residues 32–38) containing the U-turn motif (U33) and anticodon triplet (34–36). |
| **Chemical State** | Unmodified nucleotides | `[PROMPT MANDATE]` | Standard unmodified ribose-phosphate-base units (A, C, G, U). Native post-transcriptional modifications are omitted per benchmark protocol. |
| **Termini** | Standard 5'-OH / 3'-OH | `[SCIENTIFIC ASSUMPTION]` | Prevents unphysical terminal charge repulsion on closing base pair 27:43 in an excised hairpin fragment. |

---

## 4. Molecular Dynamics Simulation Specifications

All simulations will be conducted in explicit solvent under strictly controlled physiological thermodynamic conditions:

- **Force Field:** AMBER OL3 (`ff99bsc0chiOL3`). `[PROMPT MANDATE]` The OL3 parameter set specifically corrects RNA $\chi$ glycosidic torsion biases and represents the established standard for RNA hairpins.
- **Solvent Model:** Explicit water. `[PROMPT MANDATE]`
  - **Baseline Water Model:** Standard TIP3P is the default baseline because AMBER OL3 was parameterized natively in TIP3P. `[IMPLEMENTATION CHOICE]`
  - **Alternative Water Models:** Models such as OPC or TIP4P-EW remain optional only if accompanied by documented scientific justification and matched ion parameters. `[PENDING VERIFICATION]`
- **Ionic Composition:** Approximately 0.15 M potassium chloride ($\text{KCl}$): `[PROMPT MANDATE]`
  - Net neutralizing $\text{K}^+$ counterions (16 counterions for the 16 negatively charged phosphodiester linkages). `[SCIENTIFIC ASSUMPTION]`
  - Excess co-ion pairs ($\text{K}^+$ and $\text{Cl}^-$) calculated to achieve physiological $0.15\ \text{M}$ ionic strength. `[PROMPT MANDATE]`
  - **Parameter Coupling:** Ion parameters must strictly match the selected water model (e.g., `ionsjc_tip3p` for TIP3P). `[IMPLEMENTATION CHOICE]`
- **Thermodynamic Ensemble:** Isothermal-isobaric ensemble ($NPT$). `[PROMPT MANDATE]`
- **Temperature:** $T = 310.0\ \text{K}$ ($37.0\ ^\circ\text{C}$, physiological human body temperature). `[PROMPT MANDATE]`
- **Pressure:** $P = 1.0\ \text{bar}$. `[PROMPT MANDATE]`
- **Integration Timestep:** $2.0\ \text{fs}$ with constraints on all bonds involving hydrogen (SHAKE / SETTLE). `[IMPLEMENTATION CHOICE]`
- **Sampling Budget:** `[PROMPT MANDATE]`
  - 3 independent replicates per construct ($N_{\text{WT}} = 3$, $N_{\text{suppressor}} = 3$; total $N = 6$ production simulations).
  - Minimum of $20.0\ \text{ns}$ production sampling per replicate.
  - Total aggregate sampling: minimum $120.0\ \text{ns}$.
  - Independent initial velocity distributions assigned via distinct pseudo-random seeds. `[IMPLEMENTATION CHOICE]`

---

## 5. Required Trajectory Analyses

1. **Backbone RMSD (Stem vs. Loop):** `[PROMPT MANDATE]`
   - Calculated separately for the 5-bp stem (residues 27–31, 39–43) and 7-nt loop (residues 32–38).
   - Stem alignment will establish whether loop deviations represent intrinsic conformational deformation or rigid-body rotation relative to the helical axis.
2. **Per-Residue RMSF:** `[PROMPT MANDATE]`
   - Root-mean-square fluctuations of backbone atoms ($P, O5', C5', C4', C3', O3'$) computed over production trajectories to map residue-level flexibility profiles.
3. **Stem Base-Pair Hydrogen-Bond Occupancy:** `[PROMPT MANDATE]`
   - Time-resolved geometric evaluation (distance $\le 3.5\ \text{Å}$, donor–H–acceptor angle $\ge 120^\circ$) across all five Watson-Crick base pairs to verify stem stability.
4. **Canonical U-Turn Diagnostic (U33 / Phosphate 36):** `[PROMPT MANDATE]`
   - Specific hydrogen-bonding distance between U33 N3 (and/or O2') and the non-bridging phosphate oxygens of residue 36.
   - Backbone torsion angles ($\alpha_{34}, \zeta_{33}$) characterizing the sharp chain reversal between residues 33 and 34.
5. **3'-Stacking of Anticodon Bases:** `[PROMPT MANDATE]`
   - Planar stacking geometry (center-of-mass distance and normal vector angles) across pairs 34–35, 35–36, and 36–37, assessing whether the anticodon maintains the canonical stacked helical conformation.
6. **Additional Justified Structural Metric:** `[PROMPT MANDATE]`
   - One supplementary physical observable with explicit scientific justification documented prior to execution.
   - **Status:** Pending Phase 2 selection and justification. Candidates under consideration: Loop Radius of Gyration ($R_{g,\text{loop}}$), Anticodon SASA, or Stem-Loop Interhelical Tilt Angle. `[PENDING VERIFICATION]`
7. **Replicate-Aware Statistical Synthesis:** `[PROMPT MANDATE]`
   - **Sample Size Limitation:** $N=3$ independent replicates per group provide limited statistical power for traditional significance tests (e.g., two-sided Mann-Whitney $U$ test cannot reach $p < 0.05$; minimum achievable $p = 0.10$).
   - **Reporting Protocol:** Statistical synthesis will center on replicate-level effect sizes (Hedges' $g$ / Cohen's $d$, difference in replicate means/medians) with 95% bootstrap confidence intervals.
   - **Avoidance of Pseudoreplication:** Trajectory frames will not be treated as independent replicates; temporal autocorrelation will be corrected via block averaging and statistical inefficiency estimation ($\tau_{\text{int}}$).

---

## 6. Required Scientific Reporting Standards

The final scientific report must adhere to rigorous academic publishing standards:

- **Executive Conclusion:** The core research question must be answered unambiguously in the opening paragraph of the report based strictly on observed data. `[PROMPT MANDATE]`
- **Timescale Sufficiency Appraisal:** Explicit discussion of whether $3 \times 20\ \text{ns}$ ($60\ \text{ns}$ total per construct) provides sufficient phase-space sampling for RNA hairpins or whether results are subject to kinetic trapping. `[PROMPT MANDATE]`
- **Honest Characterization:** If trajectories show marginal divergence, high replicate variance, or lack of convergence, results must be characterized honestly as inconclusive rather than over-interpreted. `[PROMPT MANDATE]`
- **Translational Scope Limitations:** Explicitly articulate why isolated ASL simulations in solvent *cannot* predict cellular suppression efficiency (omission of ribosome decoding center A-site, elongation factors, codon duplex energetics, and full tRNA architecture). `[PROMPT MANDATE]`
- **Experimental Validation Pathways:** Identify specific biophysical experiments needed to validate or refute the simulation findings ($T_m$ melting profiles, NMR chemical shift mapping/NOEs, ribosome toe-printing, and cell-based dual-luciferase assays). `[PROMPT MANDATE]`
- **Omitted Modifications Impact:** Comprehensive biochemical review of omitted natural modifications at positions 34, 37, and 38, detailing their known structural roles in stabilizing the U-turn and anticodon presentation. `[PROMPT MANDATE]`
