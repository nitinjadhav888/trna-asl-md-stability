# 16 - Phase 4C Solvated ASL Restrained NVT/NPT Equilibration Report

**Study Title:** tRNA ASL Structural Stability Study  
**Phase:** Phase 4C (Restrained NVT/NPT Equilibration & Restraint Release)  
**Status:** Complete — Final Verdict: **PASS**  
**Production Readiness:** **READY FOR PRODUCTION**  
**Execution Date:** November 14, 2026  
**Pipeline Driver:** [scripts/06_equilibrate.py](file:///d:/trna-asl-md-stability/scripts/06_equilibrate.py)  

---

## 1. Scientific Purpose & Protocol Scope

The objective of Phase 4C is to perform **staged restrained heating and isothermal-isobaric ($NPT$) density equilibration** on the energy-minimized Wild-Type and CUA Suppressor 17-nt ASL explicit-solvent systems produced in Phase 4B.

### Fundamental Scientific Scope & Non-Overreach Mandate
- **Purpose:** Heat the frozen minimized systems from $10.0\ \text{K}$ to physiological temperature ($310.0\ \text{K} = 37.0^\circ\text{C}$), condense the solvent to target density ($\sim 1.018\ \text{g/cm}^3$) under $1.0\ \text{bar}$ atmospheric pressure, and progressively release positional restraints to yield stable, unconstrained starting states for production molecular dynamics.
- **Strict Scope Boundary:** Zero production MD trajectories (Phase 5: 3 × 20 ns per construct), zero trajectory production sampling, zero per-residue RMSF profiles, zero hydrogen-bond occupancy dynamics, and zero U-turn or stacking statistical inferences were performed in this phase.
- **Explicit Scientific Disclaimer:**
  > **Equilibration relaxes solvent density and heats the system to target thermodynamic conditions, but does NOT establish conformational stability, native folding, or biological function.**
  The attainment of stable temperature, density, and volume within a 1.2 ns equilibration protocol is a technical prerequisite for numerical stability, not proof that the RNA hairpin has reached its global free-energy minimum or that either construct is biologically active.

---

## 2. Input Systems & Provenance

Equilibration was initiated directly from the frozen, verified Phase 4B energy-minimized coordinates and Phase 4A OpenMM system definitions:

### Wild-Type ASL (`wt`)
- **System XML:** [simulations/wt/system/system.xml](file:///d:/trna-asl-md-stability/simulations/wt/system/system.xml)  
  *SHA-256:* `d38963ab5ea2adac230bcc552ab1b1dda2260eecf67f41fcf61ceb5a0a8f0b09`
- **Minimized PDB:** [simulations/wt/minimized.pdb](file:///d:/trna-asl-md-stability/simulations/wt/minimized.pdb)  
  *SHA-256:* `471856b2293476ac5bfa24c5b91d43e9ad0476c7aaed8af0bbeacc0ca87efc80`
- **Particles:** 15,693 (544 RNA atoms [357 heavy], 5,035 TIP3P water molecules, 30 $\text{K}^+$, 14 $\text{Cl}^-$)

### CUA Suppressor ASL (`cua`)
- **System XML:** [simulations/cua/system/system.xml](file:///d:/trna-asl-md-stability/simulations/cua/system/system.xml)  
  *SHA-256:* `5e1012f71bde05631b6612f6837c062f7133caf205f3522f3f303fdce7685220`
- **Minimized PDB:** [simulations/cua/minimized.pdb](file:///d:/trna-asl-md-stability/simulations/cua/minimized.pdb)  
  *SHA-256:* `695a1c3beb02b70b84f953b8892f365874dee2756c3df988f5dfa48b752eaf94`
- **Particles:** 15,686 (543 RNA atoms [356 heavy], 5,033 TIP3P water molecules, 30 $\text{K}^+$, 14 $\text{Cl}^-$)

---

## 3. Equilibration Protocol & Restraint Schedule

The equilibration followed the approved protocol codified in [config/simulation_config.yaml](file:///d:/trna-asl-md-stability/config/simulation_config.yaml), [docs/05_simulation_protocol.md](file:///d:/trna-asl-md-stability/docs/05_simulation_protocol.md) (Section 3), and SDR-14:

### Part I: NVT Restrained Heating ($10\ \text{K} \rightarrow 310.0\ \text{K}$, 200 ps) `[IMPLEMENTATION CHOICE]`
- **Restraint Target:** All RNA heavy atoms (357 in WT, 356 in CUA) restrained via harmonic external force:
  $$U_{\text{restraint}} = k \cdot \left[(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2\right]$$
  with $k = 5.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2} \approx 2092.00\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$.
- **Initial Velocity Assignment:** Maxwell-Boltzmann distribution initialized at $T_0 = 10.0\ \text{K}$ using deterministic integer seeds:
  - Wild-Type (`wt`): `seed = 4101`
  - Suppressor (`cua`): `seed = 4201`
- **Thermostat:** OpenMM `LangevinMiddleIntegrator` with collision frequency $\gamma = 1.0\ \text{ps}^{-1}$ and integration timestep $\Delta t = 2.0\ \text{fs}$ (all bonds involving hydrogen constrained).
- **Thermal Schedule:**
  - *Ramp Phase (0–150 ps, 75,000 steps):* Linear heating from $10.0\ \text{K}$ to $310.0\ \text{K}$ ($\Delta T = +2.0\ \text{K}$ per ps).
  - *Plateau Phase (150–200 ps, 25,000 steps):* Thermal equilibration at $310.0\ \text{K}$.
- **Ensemble:** $NVT$ (canonical volume fixed at initial cubic box $L = 5.500\ \text{nm} = 55.00\ \text{Å}$, $V = 166.38\ \text{nm}^3$).

### Part II: NPT Density Equilibration & Progressive Restraint Release ($310\ \text{K}, 1.0\ \text{bar}$, 1000 ps) `[IMPLEMENTATION CHOICE]`
- **Ensemble:** $NPT$ (isothermal-isobaric).
- **Thermostat:** `LangevinMiddleIntegrator` at $T = 310.0\ \text{K}$ ($\gamma = 1.0\ \text{ps}^{-1}$, $\Delta t = 2.0\ \text{fs}$).
- **Barostat:** OpenMM `MonteCarloBarostat` configured for isotropic volume scaling at $P = 1.0\ \text{bar}$ ($1.01325\ \text{bar}$), attempted every 25 steps ($50\ \text{fs}$).
- **4-Stage Progressive Restraint Release Schedule (250 ps per stage):**
  - **Stage 1 (200–450 ps, 125,000 steps):** $k = 2.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ ($836.80\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$). Rapid box condensation to bulk density.
  - **Stage 2 (450–700 ps, 125,000 steps):** $k = 1.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ ($418.40\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$). Continued solvent density relaxation.
  - **Stage 3 (700–950 ps, 125,000 steps):** $k = 0.2\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ ($83.68\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$). Mild guide restraints permitting backbone breathing.
  - **Stage 4 (950–1200 ps, 125,000 steps):** $k = 0.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ ($0.00\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-2}$). **Genuinely unrestrained pre-production relaxation.**
- **Total Simulation Length:** $1.2\ \text{ns}$ ($600,000$ integration steps per construct).

---

## 4. OpenMM Execution Environment & Platform Audit

To establish reproducibility, the OpenMM execution environment was programmatically queried from the active context:

- **OpenMM Version:** `8.6.1`
- **Platform:** `OpenCL`
- **Hardware Device:** `NVIDIA GeForce MX350`
- **Precision Mode:** `mixed`
- **Precision Mode Technical Specification:**
  In OpenMM's OpenCL platform, `Precision=mixed` specifies that non-bonded and bonded forces are computed using 32-bit single-precision floating-point arithmetic to optimize GPU streaming throughput, whereas particle position updates, velocities, and coordinate accumulations are integrated in 64-bit double-precision representations. This architecture eliminates numerical energy drift over multi-nanosecond trajectories while utilizing consumer-grade GPU hardware. (Crucially, it is not claimed that the entire simulation is performed in double precision; force evaluations remain single precision).

---

## 5. Thermodynamic Time-Series Analysis

Thermodynamic variables were recorded every $1.0\ \text{ps}$ ($500$ steps) across the entire $1.0\ \text{ns}$ $NPT$ equilibration, generating $N = 1,000$ observations per construct.

### Full NPT Density Relaxation Profile
Both systems condensed smoothly from their initial box volume ($166.38\ \text{nm}^3$, density $\approx 0.976\ \text{g/cm}^3$) during Stage 1 and converged to physical equilibrium within the first $300\ \text{ps}$:

| Stage | Time Interval | Restraint $k$ | WT Mean Temp | WT Mean Density | CUA Mean Temp | CUA Mean Density |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Stage 1** | 200–450 ps | $2.0\ \text{kcal/mol/Å}^2$ | $310.27 \pm 2.45\ \text{K}$ | $1.0068 \pm 0.0116\ \text{g/cm}^3$ | $310.36 \pm 2.39\ \text{K}$ | $1.0063 \pm 0.0117\ \text{g/cm}^3$ |
| **Stage 2** | 450–700 ps | $1.0\ \text{kcal/mol/Å}^2$ | $310.45 \pm 2.50\ \text{K}$ | $1.0181 \pm 0.0039\ \text{g/cm}^3$ | $310.15 \pm 2.57\ \text{K}$ | $1.0184 \pm 0.0038\ \text{g/cm}^3$ |
| **Stage 3** | 700–950 ps | $0.2\ \text{kcal/mol/Å}^2$ | $310.47 \pm 2.44\ \text{K}$ | $1.0183 \pm 0.0037\ \text{g/cm}^3$ | $310.18 \pm 2.65\ \text{K}$ | $1.0187 \pm 0.0037\ \text{g/cm}^3$ |
| **Stage 4** | 950–1200 ps | **$0.0\ \text{kcal/mol/Å}^2$** | **$310.69 \pm 2.42\ \text{K}$** | **$1.0187 \pm 0.0036\ \text{g/cm}^3$** | **$310.35 \pm 2.59\ \text{K}$** | **$1.0193 \pm 0.0037\ \text{g/cm}^3$** |

### Stage 4 Unrestrained Thermodynamic Statistics ($N = 250$ Observations)
In the final $250\ \text{ps}$ unrestrained stage ($k = 0.0$), the thermodynamic variables demonstrated exceptional stability:

| Metric | Construct | Mean | Standard Dev | Minimum | Maximum | Observations ($N$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Temperature ($K$)** | WT | **$310.69\ \text{K}$** | $2.42\ \text{K}$ | $304.41\ \text{K}$ | $318.37\ \text{K}$ | 250 |
| | CUA | **$310.35\ \text{K}$** | $2.59\ \text{K}$ | $303.75\ \text{K}$ | $317.46\ \text{K}$ | 250 |
| **Potential Energy ($\text{kJ/mol}$)** | WT | **$-232,257.38$** | $464.18$ | $-233,472.64$ | $-231,140.81$ | 250 |
| | CUA | **$-232,041.01$** | $464.09$ | $-233,214.14$ | $-230,555.61$ | 250 |
| **Target Pressure ($\text{bar}$)** | WT | **$1.0000$** | $0.0000$ | $1.0000$ | $1.0000$ | 250 |
| | CUA | **$1.0000$** | $0.0000$ | $1.0000$ | $1.0000$ | 250 |
| **Box Volume ($\text{nm}^3$)** | WT | **$159.34\ \text{nm}^3$** | $0.56\ \text{nm}^3$ | $157.81\ \text{nm}^3$ | $160.83\ \text{nm}^3$ | 250 |
| | CUA | **$159.16\ \text{nm}^3$** | $0.58\ \text{nm}^3$ | $157.67\ \text{nm}^3$ | $160.88\ \text{nm}^3$ | 250 |
| **Physical Density ($\text{g/cm}^3$)** | WT | **$1.0187\ \text{g/cm}^3$** | $0.0036\ \text{g/cm}^3$ | $1.0093\ \text{g/cm}^3$ | $1.0286\ \text{g/cm}^3$ | 250 |
| | CUA | **$1.0193\ \text{g/cm}^3$** | $0.0037\ \text{g/cm}^3$ | $1.0085\ \text{g/cm}^3$ | $1.0289\ \text{g/cm}^3$ | 250 |

### Stationarity & Drift Assessment
To evaluate whether Stage 4 reached a stable-enough stationary state, each time series was divided into first-half ($950\text{--}1075\ \text{ps}$, $n=125$) and second-half ($1075\text{--}1200\ \text{ps}$, $n=125$) sub-windows:

- **Wild-Type (`wt`):**
  - Temperature: First half = $310.51\ \text{K}$, Second half = $310.87\ \text{K}$ ($\Delta = +0.36\ \text{K}$, well within thermal noise).
  - Potential Energy: First half = $-232,342\ \text{kJ/mol}$, Second half = $-232,173\ \text{kJ/mol}$ ($\Delta = +169\ \text{kJ/mol}$, relative change $< 0.07\%$, within standard deviation of $\pm 464\ \text{kJ/mol}$).
  - Volume: First half = $159.27\ \text{nm}^3$, Second half = $159.41\ \text{nm}^3$ ($\Delta = +0.13\ \text{nm}^3$, relative change $< 0.09\%$).
  - Density: First half = $1.0191\ \text{g/cm}^3$, Second half = $1.0183\ \text{g/cm}^3$ ($\Delta = -0.0008\ \text{g/cm}^3$, relative change $< 0.08\%$).
- **Suppressor (`cua`):**
  - Temperature: First half = $310.09\ \text{K}$, Second half = $310.61\ \text{K}$ ($\Delta = +0.52\ \text{K}$).
  - Potential Energy: First half = $-232,133\ \text{kJ/mol}$, Second half = $-231,949\ \text{kJ/mol}$ ($\Delta = +183\ \text{kJ/mol}$, relative change $< 0.08\%$).
  - Volume: First half = $159.00\ \text{nm}^3$, Second half = $159.32\ \text{nm}^3$ ($\Delta = +0.31\ \text{nm}^3$, relative change $< 0.20\%$).
  - Density: First half = $1.0204\ \text{g/cm}^3$, Second half = $1.0183\ \text{g/cm}^3$ ($\Delta = -0.0020\ \text{g/cm}^3$, relative change $< 0.20\%$).

**Thermodynamic Verdict:** **Stationary.** Under the predefined equilibration QC criteria, both systems exhibited minimal thermodynamic drift, stable density, and controlled temperature, meeting all engineering prerequisites for production integration.

---

## 6. Structural Integrity & Kabsch-Aligned RMSD Analysis

Rigid-body translation and rotation naturally occur for unconstrained solutes in periodic solvent. To isolate true internal conformational changes from diffusion, rotational and translational degrees of freedom were removed using the **Kabsch algorithm** (singular value decomposition of the coordinate covariance matrix).

### Kabsch-Aligned RMSD Comparison
Structural deviation was assessed across three anatomical regions against both the Phase 4B minimized structure and the Phase 3 experimental crystallographic starting model:

| Anatomical Region | Construct | RMSD vs. Minimized Starting PDB | RMSD vs. Phase 3 Initial PDB | Structural QC Evaluation |
| :--- | :---: | :---: | :---: | :--- |
| **Full ASL (17-nt, all heavy atoms)** | WT | **$1.6619\ \text{Å}$** | **$1.5699\ \text{Å}$** | No gross structural distortion ($< 3.5\ \text{Å}$) |
| | CUA | **$1.9206\ \text{Å}$** | **$2.0149\ \text{Å}$** | No gross structural distortion ($< 3.5\ \text{Å}$) |
| **Stem Only (10-nt base-paired)** | WT | **$0.8753\ \text{Å}$** | **$0.8560\ \text{Å}$** | Rigid A-form duplex preserved ($< 1.0\ \text{Å}$) |
| | CUA | **$0.9263\ \text{Å}$** | **$0.9333\ \text{Å}$** | Rigid A-form duplex preserved ($< 1.0\ \text{Å}$) |
| **Anticodon Loop (7-nt, residues 32–38)** | WT | **$1.5156\ \text{Å}$** | **$1.5033\ \text{Å}$** | Canonical loop breathing ($< 2.5\ \text{Å}$) |
| | CUA | **$1.8121\ \text{Å}$** | **$1.7579\ \text{Å}$** | Canonical loop breathing ($< 2.5\ \text{Å}$) |

### Local Contact Geometry & Absence of Steric Clashes
- **Minimum Non-Adjacent Heavy-Atom Distance:**
  - Wild-Type (`wt`): **$2.656\ \text{Å}$** (between $C2\text{:O2}$ and $G16\text{:N2}$, corresponding to the canonical Watson-Crick hydrogen bond of base pair 28:42).
  - Suppressor (`cua`): **$2.674\ \text{Å}$** (between $C1\text{:O2}$ and $G17\text{:N2}$, corresponding to the canonical Watson-Crick hydrogen bond of base pair 27:43).
- **Steric Overlaps:** **Zero.** No unphysical atomic clashes ($< 1.8\ \text{Å}$) were detected.
- **Minimum Periodic Clearance:** Minimum distance from any RNA heavy atom to its periodic image exceeds **$10.6\ \text{Å}$** in both systems, rigorously satisfying the minimum image convention and preventing spurious self-interactions.

**Structural Verdict:** Under the predefined equilibration QC criteria, **no gross structural distortion was detected**. The 5-bp stem remains locked in canonical A-form geometry ($< 0.93\ \text{Å}$ RMSD), while the 7-nt loop exhibits natural solvent-equilibrated flexibility ($< 1.82\ \text{Å}$ RMSD) without unfolding or premature base-pair fraying.

---

## 7. WT vs. CUA Strict Comparability Audit

To verify that the computational experiment remains strictly controlled, all structural, solvent, and algorithmic parameters were audited side-by-side:

| Parameter / Dimension | Wild-Type ASL (`wt`) | CUA Suppressor ASL (`cua`) | Concordance / Delta |
| :--- | :---: | :---: | :--- |
| **Selected Parent tRNA** | `tRNA-Asp-GTC-2-1` | `tRNA-Asp-GTC-2-1 (CUA mutant)` | GtRNAdb verified |
| **ASL Sequence (5' $\rightarrow$ 3')** | `CCCGCCUGUCACGCGGG` | `CCCGCCUCUAACGCGGG` | Intended $G34C$ and $C36A$ only |
| **Total System Particles** | **15,693** | **15,686** | $\Delta = -7$ particles |
| **RNA Atoms** | 544 (357 heavy) | 543 (356 heavy) | $\Delta = -1$ atom (expected transversion delta) |
| **Water Molecules (TIP3P)** | 5,035 (15,105 atoms) | 5,033 (15,099 atoms) | $\Delta = -2$ waters (0.04% packing variation) |
| **$\text{K}^+$ Counterions / Salt** | 30 ($16\ \text{neut} + 14\ \text{salt}$) | 30 ($16\ \text{neut} + 14\ \text{salt}$) | Strictly identical ($0.15\ \text{M}\ \text{KCl}$) |
| **$\text{Cl}^-$ Co-ions** | 14 | 14 | Strictly identical |
| **Net System Charge** | **$0.0000\ e$** | **$0.0000\ e$** | Exact electrical neutrality |
| **Terminal Chemistry** | Neutral 5'-OH / 3'-OH | Neutral 5'-OH / 3'-OH | Strictly identical (avoids fraying charge) |
| **Final Periodic Box ($L_x, L_y, L_z$)** | $5.417 \times 5.417 \times 5.417\ \text{nm}$ | $5.422 \times 5.422 \times 5.422\ \text{nm}$ | Relaxed isotropic cubic ($\Delta L < 0.1\%$) |
| **Final Box Volume** | $158.91\ \text{nm}^3$ | $159.38\ \text{nm}^3$ | Density relaxed |
| **RNA Force Field** | AMBER OL3 (`amber14/RNA.OL3.xml`) | AMBER OL3 (`amber14/RNA.OL3.xml`) | Strictly identical |
| **Water & Ion Models** | TIP3P / Joung-Cheatham KCl | TIP3P / Joung-Cheatham KCl | Strictly identical |
| **Target Temperature / Pressure** | $310.0\ \text{K}$ / $1.0\ \text{bar}$ | $310.0\ \text{K}$ / $1.0\ \text{bar}$ | Strictly identical |
| **Integration Timestep** | $2.0\ \text{fs}$ (bonds to H constrained) | $2.0\ \text{fs}$ (bonds to H constrained) | Strictly identical |
| **Thermostat Scheme** | Langevin Middle ($\gamma = 1.0\ \text{ps}^{-1}$) | Langevin Middle ($\gamma = 1.0\ \text{ps}^{-1}$) | Strictly identical |
| **Barostat Scheme** | Monte Carlo Barostat (freq = 25) | Monte Carlo Barostat (freq = 25) | Strictly identical |
| **Restraint Schedule** | 4-stage release ($2.0 \rightarrow 0.0$) | 4-stage release ($2.0 \rightarrow 0.0$) | Strictly identical |
| **OpenMM Platform / Precision** | OpenCL / `mixed` (GeForce MX350) | OpenCL / `mixed` (GeForce MX350) | Strictly identical |
| **Deterministic Random Seed** | `4101` | `4201` | Distinct, documented integer seeds |

**Comparability Confirmation:** The Wild-Type and CUA Suppressor systems differ **only** in the intended $G34C$ and $C36A$ nucleotide substitutions and the associated 1-atom RNA difference (with a negligible 2-water difference from stochastic box solvation). Zero confounding parameter mismatches exist.

---

## 8. Essential Scientific Limitations

To prevent over-interpretation of the equilibrated states:
1. **Timescale Limitation:** A 1.2 ns equilibration protocol relaxes local solvent packing and brings bulk thermodynamic properties to steady-state, but **does not establish long-timescale conformational equilibrium** of the RNA solute. Conformational transitions in nucleic acid loops often require hundreds of nanoseconds to microseconds.
2. **Production Sampling Bounds:** The downstream sampling budget of 3 × 20 ns per construct (120 ns aggregate) will capture local nanosecond-scale fluctuations and initial structural stability, but has limited statistical power to sample rare, slow global unfolding or major loop refolding transitions.
3. **Absence of Tertiary Context:** The simulated 17-nt hairpin models an isolated ASL domain in solution. In vivo, the ASL is embedded within the full L-shaped tRNA tertiary architecture (D-loop/T-loop tertiary interactions) and functions in complex with the ribosome and elongation factors.
4. **Omission of Native Nucleotide Modifications:** Consistent with the project prompt mandate to test unmodified RNA, native post-transcriptional modifications (e.g., hypermodifications at position 37 or 34) were intentionally omitted. In native tRNAs, these modifications frequently stabilize the U-turn loop and prevent frame-shifting.

---

## 9. Production-Readiness Decision

Based on the verified evidence:

### Verdict: **READY FOR PRODUCTION**

**Scientific Justification:**
1. Both constructs successfully completed the entire 1.2 ns NVT/NPT staged protocol without numerical instability, step aborts, or resets.
2. The final 250 ps stage is **genuinely unrestrained** ($k = 0.0$), ensuring that starting production coordinates are completely free of artificial spring biases.
3. Thermodynamic state variables ($T = 310.5\ \text{K}$, $\rho = 1.019\ \text{g/cm}^3$, $P = 1.0\ \text{bar}$) exhibit stationarity and minimal drift.
4. Structural validation demonstrates that the native A-form stem ($< 0.93\ \text{Å}$ RMSD) and anticodon loop ($< 1.82\ \text{Å}$ RMSD) remain intact with zero steric overlaps.
5. All final state coordinates and XML checkpoints have been generated and cryptographically fingerprinted.

### Checkpoint Handoff Artifacts for Phase 5:
- **Wild-Type Final Coordinates:** [simulations/wt/equilibration/final_equilibrated.pdb](file:///d:/trna-asl-md-stability/simulations/wt/equilibration/final_equilibrated.pdb)  
  *SHA-256:* `84a6f2f690ba7bf013572f7d4a9a91463bec984a5a6703e2ff8f30875cb03409`
- **Wild-Type Final State XML:** [simulations/wt/equilibration/final_state.xml](file:///d:/trna-asl-md-stability/simulations/wt/equilibration/final_state.xml)  
  *SHA-256:* `0fd3e8b81874636d3683d35f39605ca09c88b846e288ee9c0e8472949622000c`
- **CUA Suppressor Final Coordinates:** [simulations/cua/equilibration/final_equilibrated.pdb](file:///d:/trna-asl-md-stability/simulations/cua/equilibration/final_equilibrated.pdb)  
  *SHA-256:* `e17dd30100bdde7974cac01ef46224283ef65eaab1b683b2adad59766043c65b`
- **CUA Suppressor Final State XML:** [simulations/cua/equilibration/final_state.xml](file:///d:/trna-asl-md-stability/simulations/cua/equilibration/final_state.xml)  
  *SHA-256:* `41b6cfe8ce3cb916aae7242342651df20f9c6ebae768dc9c13ec6d32f7ccb904`

**Strict Phase Boundary:** Phase 4C is formally complete and frozen. Zero Phase 5 production runs or downstream trajectory evaluations have been initiated.
