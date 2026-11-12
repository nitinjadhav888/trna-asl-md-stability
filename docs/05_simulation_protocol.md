# 05 - Molecular Dynamics Simulation Protocol

This protocol defines the step-by-step procedures for system setup, minimization, equilibration, and production MD for both the wild-type and suppressor ASL constructs.

---

## 1. System Preparation & Solvation

### 1.1 Structural Assembly
1. **Coordinate Verification:** `[SCIENTIFIC ASSUMPTION]`
   - Verify coordinate integrity of the 17-nt model (residues 27 to 43).
   - Ensure all standard RNA atoms are present ($P, O5', C5', C4', C3', O3', C1', C2', O2'$ and nucleobase atoms).
   - Assign terminal standard hydroxyl groups (5'-OH on residue 27, 3'-OH on residue 43) per SDR-09 baseline.
2. **Force Field Assignment:** `[PROMPT MANDATE]`
   - Parameterize RNA using AMBER OL3 (`ff99bsc0chiOL3`).
   - Assign explicit water parameters using TIP3P as the validated baseline force field pairing. `[IMPLEMENTATION CHOICE]`
   - **Parameter Coupling:** Monovalent ion parameters must strictly match the water model (e.g., Joung-Cheatham `ionsjc_tip3p` for TIP3P). Alternative water models require documented scientific justification and matched ion sets. `[IMPLEMENTATION CHOICE]`

### 1.2 Boundary Conditions & Box Geometry
- **Geometry:** Truncated octahedron or cubic box (selected in Phase 2). `[PENDING VERIFICATION]`
- **Padding Buffer:** Minimum of $10.0\ \text{Å}$ between any solute RNA atom and the periodic boundary edge. `[IMPLEMENTATION CHOICE]`
- **Periodic Boundary Conditions (PBC):** Applied in $x, y, z$. `[PROMPT MANDATE]`

### 1.3 Neutralization & Ion Concentration
1. **Neutralizing Counterions:** `[SCIENTIFIC ASSUMPTION]`
   - The 17-nt RNA contains 16 negative charges from the phosphodiester backbone.
   - Add exactly 16 $\text{K}^+$ counterions to achieve net electrical neutrality.
2. **Physiological Salt Addition:** `[PROMPT MANDATE]`
   - Calculate solvent volume ($V_{\text{solv}}$).
   - Add co-ions ($\text{K}^+$ and $\text{Cl}^-$) to achieve a nominal salt concentration of $0.15\ \text{M}$ ($150\ \text{mM}$):
     $$N_{\text{KCl}} = \operatorname{round}\left(0.150\ \text{mol/L} \times N_A \times V_{\text{solv}}\right)$$
   - Randomly replace non-interacting water molecules with ions, enforcing a minimum inter-ion and ion-solute exclusion distance of $5.0\ \text{Å}$. `[IMPLEMENTATION CHOICE]`

---

## 2. Energy Minimization Protocol

To eliminate steric clashes and optimize initial solvent-solute contacts, minimization proceeds in two consecutive stages:

- **Stage 1: Restrained Solute Minimization** `[IMPLEMENTATION CHOICE]`
  - Apply harmonic position restraints ($k = 5.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$) to all RNA heavy atoms.
  - Minimize solvent and ions using 2,500 steps of steepest descent followed by 2,500 steps of conjugate gradient / L-BFGS.
- **Stage 2: Full System Minimization** `[IMPLEMENTATION CHOICE]`
  - Release all restraints ($k = 0.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$).
  - Minimize the entire system for 5,000 steps until maximum force is $< 10.0\ \text{kJ}\cdot\text{mol}^{-1}\cdot\text{nm}^{-1}$ ($< 0.24\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-1}$).

---

## 3. Equilibration Protocol

Equilibration relaxes solvent density and heats the system gradually to prevent thermal shock or premature RNA denaturation:

### 3.1 NVT Restrained Heating ($0 \rightarrow 310\ \text{K}$) `[IMPLEMENTATION CHOICE]`
- **Ensemble:** $NVT$ (canonical).
- **Duration:** $200\ \text{ps}$ ($100,000$ steps at $\Delta t = 2.0\ \text{fs}$).
- **Restraints:** Harmonic restraints ($k = 5.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$) on all RNA heavy atoms.
- **Thermostat:** Langevin thermostat ($\gamma = 1.0\ \text{ps}^{-1}$).
- **Temperature Schedule:** Linearly ramp temperature from $10\ \text{K}$ to $310.0\ \text{K}$ over $150\ \text{ps}$, followed by $50\ \text{ps}$ at $310.0\ \text{K}$.
- **Initial Velocities:** Assigned from Maxwell-Boltzmann distribution at $10\ \text{K}$ using a unique random seed for each replicate.

### 3.2 NPT Density Equilibration & Restraint Release ($310\ \text{K}, 1.0\ \text{bar}$) `[IMPLEMENTATION CHOICE]`
- **Ensemble:** $NPT$ (isothermal-isobaric). `[PROMPT MANDATE]`
- **Duration:** $1.0\ \text{ns}$ ($500,000$ steps at $\Delta t = 2.0\ \text{fs}$).
- **Thermostat:** Langevin thermostat at $T = 310.0\ \text{K}$ ($\gamma = 1.0\ \text{ps}^{-1}$).
- **Barostat:** Monte Carlo barostat at $P = 1.01325\ \text{bar}$, isotropic pressure coupling, attempted every 25 steps.
- **Restraint Release Schedule:**
  - $0\text{--}250\ \text{ps}$: $k = 2.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$
  - $250\text{--}500\ \text{ps}$: $k = 1.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$
  - $500\text{--}750\ \text{ps}$: $k = 0.2\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$
  - $750\text{--}1000\ \text{ps}$: $k = 0.0\ \text{kcal}\cdot\text{mol}^{-1}\cdot\text{Å}^{-2}$ (unrestrained pre-production relaxation).

---

## 4. Unrestrained Production Molecular Dynamics

- **Ensemble:** $NPT$ ($T = 310.0\ \text{K}$, $P = 1.0\ \text{bar}$). `[PROMPT MANDATE]`
- **Production Duration:** $20.0\ \text{ns}$ per replicate ($10,000,000$ steps at $\Delta t = 2.0\ \text{fs}$). `[PROMPT MANDATE]`
- **Replicate Scheme:** `[PROMPT MANDATE]`
  - Wild-type ASL: 3 replicates (`wt_rep1`, `wt_rep2`, `wt_rep3`)
  - Suppressor ASL: 3 replicates (`suppressor_rep1`, `suppressor_rep2`, `suppressor_rep3`)
  - Total production simulations: 6
  - Total aggregate sampling: $120.0\ \text{ns}$
- **Integration Parameters:**
  - Timestep: $\Delta t = 2.0\ \text{fs}$. `[IMPLEMENTATION CHOICE]`
  - Bond Constraints: All bonds involving hydrogen constrained using SHAKE / SETTLE algorithm. `[IMPLEMENTATION CHOICE]`
  - Non-bonded Treatment: Particle Mesh Ewald (PME) for long-range electrostatics. `[IMPLEMENTATION CHOICE]`
  - Real-space Cutoff: $10.0\ \text{Å}$. `[IMPLEMENTATION CHOICE]`
  - Long-range Dispersion Correction: Applied for energy and pressure. `[IMPLEMENTATION CHOICE]`
- **Output Archiving:**
  - Trajectory frames recorded every $10.0\ \text{ps}$ ($5,000$ steps), producing $2,000$ frames per replicate ($12,000$ frames total). `[IMPLEMENTATION CHOICE]`
  - Thermodynamic state (temperature, pressure, potential energy, box volume) recorded every $10.0\ \text{ps}$. `[IMPLEMENTATION CHOICE]`
  - State checkpoints written every $1.0\ \text{ns}$ to allow restart. `[IMPLEMENTATION CHOICE]`

---

## 5. Quality Assurance & Convergence Criteria

Prior to structural analysis, every replicate trajectory must pass the following stability gates:

1. **Energy Conservation:** Total potential energy shows no systematic upward drift indicative of integration instability. `[SCIENTIFIC ASSUMPTION]`
2. **Density Plateau:** Box density converges to realistic bulk water density within the first $500\ \text{ps}$ of NPT equilibration. `[SCIENTIFIC ASSUMPTION]`
3. **Temperature & Pressure Stability:** Average temperature equals $310 \pm 1.0\ \text{K}$; average pressure equals $1.0 \pm 10.0\ \text{bar}$. `[SCIENTIFIC ASSUMPTION]`
4. **PBC Integrity:** Trajectories are centered on the RNA solute and unwrapped across periodic boundaries to prevent jump discontinuities prior to downstream analysis. `[IMPLEMENTATION CHOICE]`
