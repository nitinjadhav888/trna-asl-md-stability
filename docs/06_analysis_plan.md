# 06 - Trajectory Analysis Plan

This document formalizes the mathematical definitions, geometric criteria, algorithms, and statistical procedures for evaluating the trajectory data of the wild-type and suppressor ASLs.

---

## 1. Structural Observables & Mathematical Definitions

### 1.1 Backbone RMSD (Stem vs. Loop Dissection)
Root-Mean-Square Deviation (RMSD) evaluates conformational departure from the reference starting structure:

$$\operatorname{RMSD}(t) = \sqrt{\frac{1}{N_{\text{atoms}}} \sum_{i=1}^{N_{\text{atoms}}} \|\mathbf{r}_i(t) - \mathbf{r}_i^{\text{ref}}\|^2}$$

- **Atoms Evaluated:** Backbone heavy atoms ($P, O5', C5', C4', C3', O3'$).
- **Partitioning:**
  1. **Stem RMSD:** Evaluated over residues 27–31 and 39–43 (10 residues, 5 base pairs) after optimal superposition on the stem backbone of the reference structure.
  2. **Loop RMSD (Global):** Evaluated over loop residues 32–38 (7 residues) after alignment to the **stem** backbone. This metric captures both internal loop deformation and rigid-body hinge motion of the loop relative to the stem axis.
  3. **Loop RMSD (Internal):** Evaluated over residues 32–38 after alignment strictly to the **loop** backbone, isolating pure internal conformational reconfiguration.

---

### 1.2 Per-Residue Root-Mean-Square Fluctuation (RMSF)
RMSF measures the positional variance of each nucleotide about its time-averaged coordinate:

$$\operatorname{RMSF}_k = \sqrt{\frac{1}{T} \sum_{t=1}^T \|\mathbf{r}_k(t) - \langle\mathbf{r}_k\rangle\|^2}$$

- **Atoms Evaluated:** Average across backbone heavy atoms of nucleotide $k \in [1, 17]$.
- **Reference:** Time-averaged structure $\langle\mathbf{r}\rangle$ calculated independently for each production trajectory (after structural alignment to stem backbone).
- **Reporting:** Profile plot for positions 1–17 showing construct mean $\pm$ standard error across the 3 independent replicates.

---

### 1.3 Stem Base-Pair Hydrogen-Bond Occupancy
Evaluates the stability and intactness of the 5 canonical base pairs in the anticodon stem (27:43, 28:42, 29:41, 30:40, 31:39).

- **Geometric Criteria:**
  - Donor-Acceptor Heavy Atom Distance: $d(D, A) \le 3.5\ \text{Å}$
  - Donor–Hydrogen–Acceptor Angle: $\theta(D\text{--}H\cdots A) \ge 120.0^\circ$
- **Base Pairs Evaluated:**
  - Standard Watson-Crick G:C pairs (3 hydrogen bonds: O6–N4, N1–N3, N2–O2)
  - Standard Watson-Crick A:U pairs (2 hydrogen bonds: N6–O4, N1–N3)
- **Occupancy Fraction:**
  $$\text{Occupancy}(bp) = \frac{1}{T} \sum_{t=1}^T \mathbb{I}(\text{all canonical H-bonds intact at time } t)$$
- **Fraying Detection:** A base pair is classified as frayed in frame $t$ if zero canonical hydrogen bonds are present.

---

### 1.4 Canonical U-Turn Diagnostic (U33 / Phosphate 36)
The canonical U-turn motif (Quigley & Rich, 1976) is the defining tertiary interaction stabilizing the tRNA anticodon loop. It enables the sharp chain reversal that projects anticodon bases 34–36 outwards for ribosomal codon decoding.

```
                  Loop Residues (34-36)
                     /             \
            U33 --- (sharp turn)    A36 / Base 36
             |                       |
        H-bond: N3-H...OP1/OP2       |
             |                       |
             +-----------------------+ Phosphate 36
```

- **Diagnostic Metric 1 (Hydrogen Bond Distance):**
  - Minimum distance between U33 N3 (donor) and non-bridging phosphate oxygens of residue 36 (OP1, OP2):
    $$d_{\text{U-turn}}(t) = \min\left(\|\mathbf{r}_{\text{U33(N3)}}(t) - \mathbf{r}_{\text{P36(OP1)}}(t)\|, \|\mathbf{r}_{\text{U33(N3)}}(t) - \mathbf{r}_{\text{P36(OP2)}}(t)\|\right)$$
  - Canonical intact threshold: $d_{\text{U-turn}} \le 3.5\ \text{Å}$.
- **Diagnostic Metric 2 (Ribose O2' Interaction):**
  - Distance between U33 O2' and residue 36 acceptor (N7 in purines, or backbone).
- **Diagnostic Metric 3 (Backbone Dihedral Reversal):**
  - Monitoring the signature gauche-minus / trans conformations of torsions $\zeta_{33}$ and $\alpha_{34}$ that execute the chain reversal.

---

### 1.5 3'-Stacking Geometry of Anticodon Bases
In canonical ASL architecture, the anticodon nucleobases (34, 35, 36) form a continuous, right-handed stack with the 3'-purine at position 37.

- **Stacking Pairs Evaluated:**
  1. Residue 34 $\leftrightarrow$ Residue 35
  2. Residue 35 $\leftrightarrow$ Residue 36
  3. Residue 36 $\leftrightarrow$ Residue 37
- **Geometric Criteria:**
  - **Centroid Distance ($d_{\text{stack}}$):** Distance between the geometric centers of mass of the aromatic ring systems ($\le 4.5\ \text{Å}$).
  - **Planar Normal Tilt Angle ($\theta_{\text{tilt}}$):** Angle between normal vectors of the least-squares planes fitted to the two base rings:
    $$\theta_{\text{tilt}} = \arccos(|\mathbf{n}_A \cdot \mathbf{n}_B|)$$
  - Canonical stacked threshold: $d_{\text{stack}} \le 4.5\ \text{Å}$ and $\theta_{\text{tilt}} \le 30.0^\circ$.

---

### 1.6 Additional Justified Structural Metric (Metric 6)
- **Candidate Selected for Justification:** **Anticodon Loop Radius of Gyration ($R_{g,\text{loop}}$)**
- **Definition:**
  $$R_{g,\text{loop}}(t) = \sqrt{\frac{1}{M_{\text{loop}}} \sum_{i \in \text{loop}} m_i \|\mathbf{r}_i(t) - \mathbf{r}_{\text{COM, loop}}(t)\|^2}$$
- **Scientific Justification:** The radius of gyration measures global loop expansion or collapse. Destabilization of the U-turn or loss of 3'-stacking invariably alters the loop volume. If the suppressor anticodon disrupts loop structure, $R_{g,\text{loop}}$ will exhibit significant expansion (loop swelling/fraying) or aberrant compaction (hydrophobic collapse into loop lumen).

---

## 2. Replicate-Aware Statistical Analysis Protocol

Analyzing molecular dynamics time-series data requires rigorous handling of statistical dependencies to prevent erroneous inferences.

### 2.1 Autocorrelation Correction & Effective Sample Size
Frames recorded every 10 ps are highly correlated. To determine the number of statistically independent observations:
1. Compute the normalized autocorrelation function $C(\tau)$ for each continuous observable $X(t)$:
   $$C(\tau) = \frac{\langle (X(t) - \mu)(X(t+\tau) - \mu) \rangle}{\sigma^2}$$
2. Calculate the integrated autocorrelation time ($\tau_{\text{int}}$):
   $$\tau_{\text{int}} = 1 + 2 \sum_{\tau=1}^{\tau_{\text{max}}} C(\tau)$$
3. Compute the effective sample size for each trajectory:
   $$N_{\text{eff}} = \frac{N_{\text{total}}}{2\tau_{\text{int}}}$$

### 2.2 Replicate-Aware Inference
- Replicates represent the true independent experimental units.
- Within-replicate means ($\bar{X}_{r}$) are calculated over the production trajectory.
- The grand construct mean ($\bar{X}$) and standard error of the mean (SEM) are computed across the $N=3$ replicates:
  $$\text{SEM} = \frac{s_{\text{between-replicates}}}{\sqrt{3}}$$
- **Hypothesis Testing:** Non-parametric two-tailed permutation test and Mann-Whitney $U$ test comparing the distributions of replicate means between Wild-Type and Suppressor constructs.
- **Effect Size:** Compute Cohen's $d$ or Cliff's $\delta$ with $95\%$ bootstrap confidence intervals ($10,000$ resamples).
- **Equilibration Transient Discard:** The initial $2.0\ \text{ns}$ of production are assessed for relaxation transients and discarded from steady-state statistics if drift is detected.
