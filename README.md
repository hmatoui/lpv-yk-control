# LPV-YK Control

**Smooth switching between multiple LPV controllers via Youla–Kucera parametrization.**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reference implementation of the two control architectures introduced in:

> H. Atoui, O. Sename, V. Milanés, J.-J. Martinez,
> *"Advanced Design of Multiple LPV Controllers with Youla–Kucera Parametrization: Real-world Validation on Autonomous Vehicles."*

Both were experimentally validated on a robotized Renault ZOE at the Satory test track (France).

---

## What problem this solves

A single LPV controller synthesized over a wide scheduling region is conservative, for two separate reasons:

- **Lyapunov conservatism** — one Lyapunov function must certify every operating point at once.
- **Over-bounding** — in the polytopic approach the operating region is replaced by a convex hull. When the scheduling parameters are coupled (the canonical case being $\rho_1 = v_x$ and $\rho_2 = 1/v_x$ in the lateral bicycle model), that hull contains points the plant can never reach, and the design is forced to stabilize physically impossible models.

Partitioning the region and switching between local controllers fixes both — but classical switched-LPV schemes charge a high price for it: all local controllers must be **co-designed** in one coupled LMI problem, stability holds only for **restricted switching signals** (dwell-time, hysteresis), and smooth transitions require **extra bumpless-transfer machinery**.

This repository implements two architectures that remove all three costs by routing the switching through the Youla–Kucera parameter instead of the controller itself.

| | Classical switched LPV | LPV–YK (this repo) |
|---|---|---|
| Local controller design | Co-designed under shared LMIs | **Independent**, any method (PID, H∞, LQG) |
| Admissible switching | Dwell-time / hysteresis constrained | **Any bounded signal**, continuous or not |
| Transition smoothness | Needs bumpless-transfer add-ons | **Structural** — no controller state reset |
| Adding a controller | Re-certify the whole design | **Plug & play** — add one $Q_i$ |
| LMI problem size | Grows with number of subregions | **Independent** of it |

---

## The core idea

Given a plant $G$ and one stabilizing nominal controller $K_0$, every stabilizing controller is

$$K(Q) = (U_0 + M_0 Q)\,(V_0 + N_0 Q)^{-1}$$

for some **stable** $Q$ — and the closed loop is *affine* in $Q$. So instead of certifying a switched controller, you certify a switched $Q$, where stability is structural rather than the outcome of an optimization. $Q = 0$ gives back $K_0$; each local controller $K_i$ corresponds to one $Q_i$.

After a state transformation the closed-loop state matrix becomes block upper-triangular, with the switching signal $\sigma$ appearing **only in bounded off-diagonal blocks**. Since bounded off-diagonal terms cannot destabilize a block-triangular system, switching cannot destabilize the loop — for any bounded $\sigma$.

A full derivation, including both stability proofs, is in [`docs/theory.md`](docs/theory.md).

---

## The two architectures

### Grid-based LPV–YK (`example01`)

For general LPV plants, **including parameter-varying input/output matrices**. Local LPV controllers are designed per subregion; the $N$ Youla parameters run as a parallel bank whose outputs are selected by $\sigma$:

$$Q(\rho,\sigma) = \sum_{i=1}^{N}\sigma_i(\rho)\,Q_i(\rho)$$

Because $A_q = \mathrm{diag}(A_{q,1},\dots,A_{q,N})$ and $B_q$ are **independent of $\sigma$**, every internal state keeps evolving continuously across a switch. There is no state re-initialization, which is exactly why transitions are bumpless without any added machinery.

Guarantees exponential stability with $\gamma_\infty = \max_i \gamma_{\infty,i}$ — you inherit the worst of your *local* designs, not a global compromise.

### Partitioned polytopic LPV–YK (`example02`)

For strictly proper plants with constant $B_2, C_2, D_{12}, D_{21}$. LTI controllers are designed at polytope vertices, YK-interpolated inside each subset, and switched between subsets. Requires **no rate bounds** on $\dot\rho$.

Two ideas carry this design:

1. **Trajectory-aligned partitioning.** A chain of small polytopes hugging the physical curve $\rho_2 = 1/\rho_1$ covers far less fictitious territory than one hull.
2. **Shared-face controllers (A.3.2).** Adjacent subsets share a vertex; the controller there is designed **once and reused**, so $K_{i,3} \equiv K_{i+1,2}$. Switching across that face is then not merely stable — it is a *no-op*. In the published experiments the scheduling signal chattered between two subsets at $t \in [28.3, 28.5]$ s with **no effect** on steering or yaw rate.

This repo enforces that identity structurally: vertex controllers are keyed on their coordinates and designed once, and `example02` asserts $K_{p,3} = K_{p+1,2}$.

### Choosing between them

| | Grid-based | Partitioned polytopic |
|---|---|---|
| Plant class | General LPV, varying $B_2, C_2$ | Strictly proper, constant $B_2, C_2, D_{12}, D_{21}$ |
| Local controllers | LPV, per subregion | LTI, per vertex |
| Stability certificate | Parameter-dependent $X(\rho)$ | Constant $X$ (quadratic) |
| Rate bounds $\dot\rho$ | **Required** | **Not required** |
| Offline cost | Higher, scales with grid | Lighter, scales with $2^{n_p}$ |
| Boundary behaviour | Bumpless | **Exactly identical controller** |

Use the grid-based scheme for parameter-varying input/output matrices or lower conservatism; use the partitioned polytopic scheme when over-bounding dominates or scheduling rates are fast or unknown.

---

## Repository layout

```
LPV-YK-Control/
├── setupPath.m                 Add src/ to the MATLAB path
├── checkDependencies.m         Report which toolboxes are available
├── src/
│   ├── model/                  Vehicle, actuator, bicycle model (grid + polytopic)
│   ├── weights/                H-infinity weighting filters We, Wu
│   ├── synthesis/              Controller design wrappers
│   ├── lmi/                    LMI solvers (YALMIP)
│   ├── youla/                  Coprime factorization, Youla parameter, YK controller
│   └── utils/                  Partitioning, interpolation, order padding, checks
├── examples/
│   ├── example01_gridLpvYk.m
│   └── example02_partitionedPolytopicLpvYk.m
├── tests/                      MATLAB unit tests
└── docs/theory.md              Full derivation and proof walkthrough
```

### The functions that matter

| Function | Role | Paper |
|---|---|---|
| `coprimeFactorization` / `coprimeFactorizationLpv` | Doubly-coprime factorization → the fixed block $J(\rho)$ | eq. (11), (33) |
| `youlaParameter` | $Q_i$ encoding a local controller w.r.t. the nominal | eq. (10), (32) |
| `ykController` | $(U+MQ)(V+NQ)^{-1}$, both left and right forms | eq. (8), (31) |
| `lmiStateFeedbackGrid` | Rate-bounded pLMIs → $F_g$, $F_{k,0}$ | eq. (6)–(7) |
| `lmiHinfStateFeedbackPolytope` | Vertex LMIs → $F_g$, $F_{k,0}$ | eq. (28)–(29) |
| `lmiHinfPolytope` | Polytopic LPV H∞ output-feedback synthesis | (A.3.1) |
| `polytopeVertices` | Triangular hull, drops the infeasible corner | Fig. 6 |
| `equalizeOrder` | Pad $Q_i$ to a common order for the block-diagonal bank | eq. (9) |

---

## Requirements

| Component | Grid-based | Partitioned polytopic |
|---|:---:|:---:|
| MATLAB R2019b or newer | ✔ | ✔ |
| [Control System Toolbox](https://www.mathworks.com/products/control.html) | ✔ | ✔ |
| [Robust Control Toolbox](https://www.mathworks.com/products/robust.html) | ✔ | ✔ |
| [LPVTools](https://www.mathworks.com/matlabcentral/fileexchange/59563-lpvtools) | ✔ | — |
| [YALMIP](https://yalmip.github.io/) | optional¹ | ✔ |
| [SeDuMi](https://sedumi.ie.lehigh.edu/) or [SDPT3](https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/) | optional¹ | ✔ |

¹ The grid pipeline can take $F_g$, $F_{k,0}$ either from LPVTools' `lpvsfsyn` or from the bundled `lmiStateFeedbackGrid`, which needs YALMIP.

Run `checkDependencies` to see what is installed.

---

## Quick start

```matlab
setupPath;
checkDependencies;

% Grid-based LPV-YK over [5, 30] m/s, 5 subsets
run examples/example01_gridLpvYk.m

% Partitioned polytopic LPV-YK on the trajectory-aligned triangular hull
run examples/example02_partitionedPolytopicLpvYk.m

% Unit tests (fast, needs only Control System Toolbox)
runtests('tests')
```

### Minimal use of the Youla layer

```matlab
setupPath;

% Your plant and any two stabilizing controllers
G  = ss(...);        % plant
K0 = ss(...);        % nominal / central controller, globally valid
K1 = ss(...);        % local controller, designed independently

% Factorize once, relative to the nominal controller
cf = coprimeFactorization(G, K0, Fg, Fk0);     % pass [] for pole-placement fallback

% Encode the local controller as a Youla parameter
Q1 = youlaParameter(G, K0, K1, cf.Fg, cf.Fk0);

% Realize it - and confirm the recovery is exact, not approximate
[K, Kdual] = ykController(cf, Q1);
verifyStepEquality({K1, K, Kdual}, {'K1', 'right form', 'left form'});
```

Adding a second local controller later is one more call to `youlaParameter`. Nothing already computed changes — that is the plug-and-play property.

---

## Reproducing the published results

`example01` prints the γ-performance of each local LPV controller (Table I in the paper) and `example02` prints the vertex LTI controllers' γ values (Table II), alongside the nominal controllers' levels.

Published values, for comparison:

| Design | Published γ |
|---|---|
| Nominal grid LPV controller $K_0$ over $[5,30]$ | 1.55 |
| Single grid LPV controller, full region, local weights | 1.42 |
| Local grid controllers $K_1 \dots K_5$ | 1.1543, 1.1157, 1.2150, 1.2112, 1.2175 |
| **Worst local level** $\max_i \gamma_{\infty,i}$ | **1.2175** |
| Nominal polytopic $K_0$ over $\mathcal{P}_0$ | 1.60 |
| LTI vertex controllers $K_{ij}$ | ≈ 0.97 – 1.12 |

The grid comparison — **1.2175 versus 1.42** for the same weights — is the quantitative payoff: roughly a 14% better guaranteed $\mathcal{L}_2$ level, obtained purely by partitioning.

`example02` additionally sweeps the closed-loop H∞ norm against speed for the partitioned design versus a gain-scheduled LPV–YK design over the full hull. Both use identical YK machinery, so the gap isolates the benefit of the trajectory-aligned partition.

> **Note.** Experimental data, the Simulink/dSPACE real-time models, and the reinforcement-learning extensions of the original research code are deliberately **not** included here. This repository is the control-design pipeline only.

---

## Notes on numerical behaviour

- **Order growth.** The grid-based controller has order $n_x + n_{k,0} + \sum_i n_{q,i}$, and all $N$ Youla parameters run continuously. This is the acknowledged cost of the architecture. `equalizeOrder` pads the shorter $Q_i$ with fast, decoupled modes so the block-diagonal bank has a fixed state dimension, which a real-time target requires.
- **Design freedom in $F_g$, $F_{k,0}$.** Closed-loop stability holds for *any* admissible pair, but the choice shapes the dynamics of $Q$ and hence transient behaviour during switching. Optimizing it is named as open work in the paper's conclusion. `lmiStateFeedbackRegion` constrains the factorization poles to a well-conditioned LMI region as a practical default.
- **Performance level at intermediate $\sigma$.** The bound $\gamma_\infty = \max_i \gamma_{\infty,i}$ rests on YK performance recovery, which is exact when $\sigma$ sits at a vertex of $\{0,1\}^N$. Stability holds unconditionally for any bounded $\sigma$; the γ bound is cleanest for selector-type switching, which is what the examples implement.

---

## Citation

If you use this code, please cite the paper:

```bibtex
@article{atoui_lpvyk,
  title   = {Advanced Design of Multiple {LPV} Controllers with
             {Youla}--{Kucera} Parametrization: Real-world Validation
             on Autonomous Vehicles},
  author  = {Atoui, Hussam and Sename, Olivier and Milan{\'e}s, Vicente
             and Martinez, John-Jairo},
  note    = {See CITATION.cff for the current publication details}
}
```

Related work by the authors:

- H. Atoui, O. Sename, V. Milanés, J.-J. Martinez, *"Interpolation of multi-LPV control systems based on Youla–Kucera parameterization,"* **Automatica**, 134:109963, 2021.
- H. Atoui, O. Sename, V. Milanés, J.-J. Martinez, *"Toward switching/interpolating LPV control: A review,"* **Annual Reviews in Control**, 54:49–67, 2022.
- H. Atoui, O. Sename, V. Milanés, J.-J. Martinez, *"LPV-Based Autonomous Vehicle Lateral Controllers: A Comparative Analysis,"* **IEEE T-ITS**, 2021.

---

## Acknowledgements

The LMI synthesis routines build on code by **Charles Poussot-Vassal** (`lmiHinfPolytope`), **Olivier Sename** (`lmiHinfStateFeedbackPolytope`) and **Kazusa Yamada** (`lmiStateFeedbackGrid`), developed at GIPSA-lab, Grenoble. Original authorship is preserved in each file header.

Experimental validation was supported by the RENAULT research department.

## License

[MIT](LICENSE) © Hussam Atoui
