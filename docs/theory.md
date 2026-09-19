# Advanced Design of Multiple LPV Controllers with Youla–Kucera Parametrization

**Theoretical background for this repository: a detailed walkthrough of the paper, its two contributions, and how each equation maps onto the code in `src/`.**

Atoui, Sename, Milanés, Martinez — *Advanced Design of Multiple LPV Controllers with Youla–Kucera Parametrization: Real-world Validation on Autonomous Vehicles.*

---

## Table of contents

1. [The problem the paper attacks](#1-the-problem-the-paper-attacks)
2. [Background: why single-Lyapunov LPV design runs out of room](#2-background-why-single-lyapunov-lpv-design-runs-out-of-room)
3. [The Youla–Kucera idea in one page](#3-the-youlakucera-idea-in-one-page)
4. [Notation](#4-notation)
5. [Contribution I — Grid-based LPV–YK switching control](#5-contribution-i--grid-based-lpvyk-switching-control)
6. [Contribution II — Partitioned polytopic-based LPV–YK control](#6-contribution-ii--partitioned-polytopic-based-lpvyk-control)
7. [Lemma A.1 — the triangular-stability workhorse](#7-lemma-a1--the-triangular-stability-workhorse)
8. [Side-by-side comparison of the two architectures](#8-side-by-side-comparison-of-the-two-architectures)
9. [Application: lateral control of an autonomous vehicle](#9-application-lateral-control-of-an-autonomous-vehicle)
10. [Experimental results](#10-experimental-results)
11. [Mapping theory to the code in this repository](#11-mapping-theory-to-the-code-in-this-repository)
12. [Critical reading: strengths, limits, editorial notes](#12-critical-reading-strengths-limits-editorial-notes)

---

## 1. The problem the paper attacks

An LPV controller synthesized over a **wide** scheduling region is conservative. Two distinct mechanisms cause this:

- **Lyapunov conservatism.** One Lyapunov function (constant, or parameter-dependent with a fixed basis) must certify every operating point at once. Over a large region such a function may be very poor, or may not exist.
- **Over-bounding of the parameter set.** In the polytopic approach the operating region is replaced by a convex hull. When parameters are *coupled* — the canonical automotive case being $\rho_1 = v_x$ and $\rho_2 = 1/v_x$ in the lateral bicycle model — the hull contains vertices that the physical system never visits, and the design is forced to stabilize plants that do not exist. Some of those fictitious vertex models may even be unstable.

The standard cure is **switched LPV control**: partition the region, design a controller per subregion, and switch. The pioneering work here is Lu & Wu (ref. [26]), extended for smoothness by Hanifzadegan & Nagamune [16] and He et al. [17], [39].

The **price** of those methods is what this paper targets:

1. All local controllers must be **co-designed** in one coupled LMI/PLMI problem, because the switching stability certificate couples them. Adding a sixth controller to a five-controller design means redoing the whole synthesis.
2. Stability holds only for **restricted switching signals** — hysteresis switching, average dwell-time, minimum dwell-time.
3. Smoothness at switching instants requires **extra machinery**: interpolation functions, higher-order differential control signals, iterative descent over three coupled decision variables, or explicit bumpless-transfer mechanisms.

The paper's claim is that Youla–Kucera parametrization removes all three costs simultaneously, and it backs this with real-car experiments.

**The four advertised advantages:**

| # | Advantage | What it replaces |
|---|---|---|
| 1 | Local LPV controllers designed **separately**; LMI conditions decoupled | Coupled co-design PLMIs |
| 2 | **Arbitrary** bounded switching signals — continuous or discontinuous | Dwell-time / hysteresis constraints |
| 3 | Smooth (bumpless) control and state responses at switching | Explicit bumpless-transfer add-ons |
| 4 | **Plug & play**: add or remove a local controller without re-certifying | Full re-synthesis |

---

## 2. Background: why single-Lyapunov LPV design runs out of room

### 2.1 Grid-based approach

Solve parameter-dependent LMIs (pLMIs) on a grid of $\rho$ using a single smooth parameter-dependent Lyapunov function $X(\rho)$. Handles non-affine parameter dependence naturally (you never need affinity — you just evaluate at grid points), and exploits rate bounds $\underline{\nu}_k \le \dot\rho_k \le \overline{\nu}_k$. Its limitation is the single $X(\rho)$ over a large region.

Lu & Wu's fix: multiple parameter-dependent Lyapunov functions $X_i(\rho)$, one per subregion, plus junction conditions at the switching surfaces ensuring the Lyapunov function does not jump upward when switching. Those junction conditions are exactly what couples the designs.

### 2.2 Polytopic approach

Write the plant as a convex combination of vertex LTI models, solve LMIs at the vertices, use a **constant** Lyapunov matrix (quadratic stability). The most popular approach — and the most conservative, for two compounding reasons: quadratic stability, plus over-bounding.

Attempts in the literature to reduce over-bounding: trajectory-aligned convex polyhedrons (leads to vertex explosion), Scheduling Dimension Reduction via PCA (Kwiatkowski & Werner [23]), DNN-based SDR (Koelewijn & Tóth [22]), overlapped-subset switching (Zhao & Nagamune [44], [45]).

### 2.3 Prior YK-based interpolation work

- Rasmussen & Chang [33]: YK configuration to improve a polytopic LPV control.
- Bianchi & Peña [10]: YK-based gain-scheduling of **LTI** controllers designed at the vertices of **one** polytope. Quadratic stability and performance guaranteed at interpolation points. **This paper's Lemma III.1 is theirs.**
- Blanchini et al. [11]: fixed pole-assignment by YK interpolation.
- Atoui et al. [6] (the authors' own *Automatica* 2021 paper): interpolation between **two** LPV controllers, each designed over the **full** polytopic region.

**The gap this paper fills:** all of the above interpolate controllers defined over a *single* region (LTI controllers at the vertices of one polytope, or LPV controllers over the whole polytope). Switching among **multiple LPV controllers designed over distinct parameter subregions** had not been addressed. That is the novelty claim, and it is what makes the partitioning — and hence the reduction in over-bounding — possible.

---

## 3. The Youla–Kucera idea in one page

Given a plant $G$ and one stabilizing controller $K_0$ (the *nominal* or *central* controller), the set of **all** controllers that stabilize $G$ is parametrized by a single stable transfer function $Q$:

$$
K(Q) = (U_0 + M_0 Q)(V_0 + N_0 Q)^{-1} = (\tilde U_0 + Q\tilde M_0)(\tilde V_0 + Q \tilde N_0)^{-1}
$$

where $G = N_0 M_0^{-1} = \tilde M_0^{-1}\tilde N_0$ and $K_0 = U_0 V_0^{-1} = \tilde V_0^{-1}\tilde U_0$ are right/left coprime factorizations over $\mathcal{RH}_\infty$.

Three consequences drive the whole paper:

1. **Stability becomes a property of $Q$ alone.** $K(Q)$ stabilizes $G$ **if and only if** $Q$ is stable. The closed-loop map is *affine* in $Q$. So instead of certifying a switched controller, you certify a switched $Q$ — and $Q$ lives in a space where stability is structurally easy.
2. **$Q = 0 \Leftrightarrow K = K_0$.** The nominal controller is the origin of the parameter space. Every other controller is a stable "perturbation" away from it.
3. **Performance recovery.** A parametrized controller $\tilde K(Q)$ recovers the closed-loop performance of the actual controller $K$ it encodes (Theorem 3.9 in Navas [28]). So $\mathcal{F}_l(G, \tilde K) \equiv \mathcal{F}_l(G, K)$ — realizing $K_i$ through the YK structure costs nothing in performance.

The paper implements $K(Q)$ as a **lower LFT** $\tilde K = \mathcal{F}_l(J, Q)$, where $J$ is a fixed observer-based generalized controller built from the nominal $K_0$ plus two state-feedback gains, and $Q$ is the switched/interpolated part. Switching then touches **only** $Q$ — $J$ stays frozen. This is the structural reason the local designs decouple.

The two contributions are two different ways of constructing $Q$:

- **Grid-based**: $Q(\rho,\sigma) = \sum_{i=1}^N \sigma_i(\rho)\, Q_i(\rho)$ — a *parallel bank* of $N$ LPV YK parameters, output-selected by the switching vector $\sigma$.
- **Partitioned polytopic**: $\bar Q_i(\rho) = \sum_j \alpha_{ij}(\rho)\bar Q_{ij}$ — a *convex interpolation* of vertex YK parameters inside each subset, with switching between subsets.

---

## 4. Notation

| Symbol | Meaning |
|---|---|
| $\mathcal{P}$ | full scheduling parameter region; $\mathcal{P} = \bigcup_i \mathcal{P}_i$ |
| $\mathcal{P}_i$, $i\in\mathbb{I}_N=\{1..N\}$ | closed parameter subsets (subregions) |
| $\mathcal{P}_0$ | full region (grid case) / convex hull containing all trajectories (polytopic case) |
| $K_0(\rho)$ | nominal (central) LPV controller, valid over the whole region |
| $K_i(\rho)$ | local LPV controller for subregion $\mathcal{P}_i$ |
| $K_{ij}$ | local **LTI** controller at vertex $j$ of subset $i$ (polytopic case) |
| $w_{ij}$ | $j$-th vertex of polytope $\mathcal{P}_i$, $j \in \mathbb{I}[1, 2^{n_p}]$ |
| $\alpha_{ij}(\rho)$ | polytopic scheduling coefficients, $\sum_j \alpha_{ij}=1$, $\alpha_{ij}\ge0$ |
| $\sigma = [\sigma_1 \dots \sigma_N]$ | switching vector, $\sigma_i \in [0,1]$ |
| $Q_i(\rho)$ | YK parameter encoding $K_i$ relative to $K_0$ |
| $F_g(\rho)$ | plant state-feedback gain used to build the coprime factors |
| $F_{k,0}(\rho)$ | nominal-controller state-feedback gain |
| $J(\rho)$ | fixed part of the LFT controller realization |
| $\gamma_{\infty,i}$ | $\mathcal{L}_2$-induced performance level of local loop $i$ |
| $\{\underline\nu_k, \overline\nu_k\}$ | the LMI must hold at **both** rate-bound extremes (all vertices of the rate box) |

Subscript $i$ = subregion index; second subscript $j$ = polytope vertex index. E.g. $A_{K,ij}$ = state matrix of the $i$-th local controller at the $j$-th vertex.

---

## 5. Contribution I — Grid-based LPV–YK switching control

### 5.1 Plant and standing assumptions

MIMO LPV plant, eq. (1):

$$
G(\rho):\quad
\begin{cases}
\dot x = A(\rho)x + B_1(\rho)w + B_2(\rho)u\\
z = C_1(\rho)x + D_{11}(\rho)w + D_{12}(\rho)u\\
y = C_2(\rho)x + D_{21}(\rho)w + D_{22}(\rho)u
\end{cases}
$$

with $w = [r\ \ n\ \ d]^T$ (reference, noise, input disturbance). Standing assumptions:

- $(A(\rho), B_2(\rho), C_2(\rho))$ parameter-dependent stabilizable and detectable $\forall\rho\in\mathcal{P}$;
- $[B_2^T(\rho)\ \ D_{12}^T(\rho)]$ and $[C_2(\rho)\ \ D_{21}(\rho)]$ full row rank;
- $D_{22}(\rho) = 0$ (strictly proper input-to-measured-output);
- $\rho$ in a compact set with bounded rates $\underline\nu_k \le \dot\rho_k \le \overline\nu_k$.

**Note the generality here:** $B_2$, $C_2$, $D_{12}$, $D_{21}$ are all allowed to be parameter-varying. This is the grid-based approach's structural advantage and the reason it is presented first.

**Partition.** $\mathcal{P}$ is covered by finitely many closed subsets $\{\mathcal{P}_i\}$ separated by switching surfaces. Two subsets are *adjacent* if their closures intersect on a nonempty boundary. **Overlaps are not required** and are allowed only on boundaries — a notable simplification versus overlapped-subset switching schemes.

### 5.2 Assumption (A.2.1) — the decoupling assumption

> Over each subset $\mathcal{P}_i$, $i \in \mathbb{I}_N \cup \{0\}$, there exists an LPV controller $K_i(\rho)$, **pre-designed separately**, that exponentially stabilizes $G(\rho)$ over $\mathcal{P}_i$ and achieves a suitable performance there.

$$
K_i(\rho): \begin{bmatrix} A_{K,i}(\rho,\dot\rho) & B_{k,i}(\rho)\\ C_{k,i}(\rho) & D_{k,i}(\rho)\end{bmatrix}, \quad i \in \{0\}\cup\mathbb{I}_N
$$

This is the entire interface to the local designs. Nothing about $K_i$ beyond "it stabilizes $G$ on $\mathcal{P}_i$ with level $\gamma_{\infty,i}$" is used. Hence advantage 4 in the table above: the local controllers can even come from **different synthesis paradigms** — PID, $H_\infty$, LQG — as long as each is stabilizing on its own patch.

$K_0$ plays a distinguished role: it must be valid over the **whole** region $\mathcal{P}_0$. It may be conservative — that is fine and expected. It is the "safety net" whose performance the local $Q_i$'s improve upon.

Each local closed loop $CL_i(\rho)$ satisfies the bounded real lemma with level $\gamma_{\infty,i}$ via **Lemma II.1** (eq. 4), the standard rate-bounded parameter-dependent BRL:

$$
\begin{bmatrix}
A_{cl,i}^T X_{cl,i} + X_{cl,i}A_{cl,i} + \sum_k \{\underline\nu_k,\overline\nu_k\}\tfrac{\partial X_{cl,i}}{\partial\rho_k}
& X_{cl,i}B_{cl,i} & C_{cl,i}^T\\
B_{cl,i}^T X_{cl,i} & -\gamma_{\infty,i}I_{n_w} & D_{cl,i}^T\\
C_{cl,i} & D_{cl,i} & -\gamma_{\infty,i}I_{n_z}
\end{bmatrix} < 0
$$

### 5.3 The switching logic (A.2.2)

$$
\sigma_i(\rho) = 0\ \forall i \;\Longrightarrow\; \tilde K(\rho,\sigma) \equiv K_0(\rho)
$$
$$
\rho \in \mathcal{P}_l \;\Longrightarrow\; \sigma_l(\rho) = 1,\ \ \sigma_i(\rho) = 0\ \forall i \ne l \;\Longrightarrow\; \tilde K(\rho,\sigma) \equiv K_l(\rho)
$$

So $\sigma$ is a **selector**. In the vehicle application it is a hard $\{0,1\}$ switch taken when $\rho(t)$ hits a boundary. The theorem, however, is proved for **any bounded $\sigma_i \in [0,1]$** — continuous, discontinuous, chattering, hysteretic. That generality is the point: stability does not depend on how $\sigma$ behaves.

**Remark II.1** notes that with multiple parameters the transition is always to a *unique* adjacent subregion — each $\rho$ belongs to exactly one subset, hence to exactly one controller — so the structure does not change with $\dim\rho$.

### 5.4 Theorem II.1 — the synthesis conditions

> Under (A.2.1), the parametrized LPV–YK controller $\tilde K(\rho,\sigma)$ of (5)–(8) exponentially stabilizes $G(\rho)$ for **any** continuous/discontinuous bounded $\sigma_i$, if there exist symmetric positive-definite parameter-dependent $X_g(\rho)$, $X_{k,0}(\rho)$ and matrices $V(\rho)$, $W(\rho)$ such that $\forall\rho\in\mathcal{P}$:

$$
A(\rho)X_g + X_g A^T(\rho) - \sum_{j=1}^p \{\underline\nu_j,\overline\nu_j\}\frac{\partial X_g}{\partial\rho_j} + B_2(\rho)V(\rho) + V^T(\rho)B_2^T(\rho) < 0 \tag{6}
$$

$$
A_{K,0}(\rho)X_{k,0} + X_{k,0}A_{K,0}^T(\rho) - \sum_{j=1}^p \{\underline\nu_j,\overline\nu_j\}\frac{\partial X_{k,0}}{\partial\rho_j} + B_{k,0}(\rho)W(\rho) + W^T(\rho)B_{k,0}^T(\rho) < 0 \tag{7}
$$

with $F_g(\rho) = V(\rho)X_g^{-1}(\rho)$ and $F_{k,0}(\rho) = W(\rho)X_{k,0}^{-1}(\rho)$.

**What these two LMIs actually are.** Substituting $V = F_g X_g$ turns (6) into

$$
(A + B_2F_g)X_g + X_g(A+B_2F_g)^T - \sum_j\{\underline\nu_j,\overline\nu_j\}\frac{\partial X_g}{\partial \rho_j} < 0,
$$

i.e. the textbook rate-bounded LPV **state-feedback stabilizability** condition in the "$X$-form". (6) says: the *plant* is stabilizable by parameter-dependent state feedback. (7) says the same for the *nominal controller's* dynamics $A_{K,0}$ — this is what **Remark II.2** flags, and it is the usual requirement for constructing a coprime factorization of $K_0$ (cf. Bianchi & Peña [10], Hespanha & Morse [18]).

The sign convention flips between (6) and the Lyapunov inequality (14a) because $Y_g = X_g^{-1}$ and $\partial(X^{-1}) = -X^{-1}(\partial X)X^{-1}$ — consistent, but worth noticing when coding it.

**Crucially: neither (6) nor (7) mentions $K_i$, $\sigma$, or $N$.** The size of the synthesis problem is *independent of the number of subregions*. That is the formal statement of advantages 1 and 4.

### 5.5 The controller realization (eq. 8)

$\tilde K(\rho,\sigma)$ has three block states: a plant-state observer-like copy, the nominal controller state $x_{K,0}$, and the stacked YK-parameter state $x_q$.

$$
\tilde A_K(\rho,\sigma) = \begin{bmatrix}
A + B_2F_g - B_2D_qC_2 & -B_2D_qF_{k,0} & B_2C_q\\
-B_{k,0}C_2 & A_{K,0} & 0\\
-B_qC_2 & -B_qF_{k,0} & A_q
\end{bmatrix},\quad
\tilde B_k = \begin{bmatrix}B_2D_q\\ B_{k,0}\\ B_q\end{bmatrix}
$$

$$
\tilde C_k = \begin{bmatrix} F_g - (D_{k,0}+D_q)C_2 & C_{k,0}-D_qF_{k,0} & C_q\end{bmatrix},\qquad
\tilde D_k = D_{k,0} + D_q
$$

with the **switched YK bank** (eq. 9):

$$
A_q = \mathrm{diag}(A_{q,1},\dots,A_{q,N}),\quad
B_q = \begin{bmatrix}B_{q,1}^T \cdots B_{q,N}^T\end{bmatrix}^T
$$
$$
C_q(\rho,\sigma) = \begin{bmatrix}\sigma_1 C_{q,1} & \cdots & \sigma_N C_{q,N}\end{bmatrix},\qquad
D_q(\rho,\sigma) = \sum_{i=1}^N \sigma_i D_{q,i}
$$

**Read the structure carefully — this is the mechanism.** All $N$ YK parameters run **in parallel, always, all the time**. $A_q$ and $B_q$ do not depend on $\sigma$: every $Q_i$ is continuously driven by the same input, and all $N$ internal states $x_{q,i}$ evolve continuously regardless of the switching. $\sigma$ appears **only in the output map** $C_q$ and the feedthrough $D_q$.

This is why switching is bumpless without any bumpless-transfer machinery: there is **no controller state re-initialization at a switch**. The state vector is continuous across switching instants by construction; only the weighting of already-running, already-stable dynamics changes. Contrast with classical switched LPV, where the incoming controller's state must be initialized somehow (hence state-reset schemes [9], bumpless transfer [15], [21], [35]).

The cost is order: $n_x + n_{k,0} + \sum_i n_{q,i}$. For $N=5$ this is a large controller. The paper concedes this ("even though they result in high order state-space controllers").

### 5.6 The YK parameter $Q_i(\rho)$ (eq. 10)

$$
A_{q,i} = \begin{bmatrix}
A + B_2D_{k,i}C_2 & B_2C_{k,i} & B_2[D_{k,i}-D_{k,0}]F_{k,0} - B_2C_{k,0}\\
B_{k,i}C_2 & A_{K,i} & B_{k,i}F_{k,0}\\
0 & 0 & A_{K,0}+B_{k,0}F_{k,0}
\end{bmatrix}
$$
$$
B_{q,i} = \begin{bmatrix}B_2[D_{k,i}-D_{k,0}]\\ B_{k,i}\\ B_{k,0}\end{bmatrix},\qquad
D_{q,i} = D_{k,i} - D_{k,0}
$$
$$
C_{q,i} = \begin{bmatrix}D_{k,i}C_2 - F_g & C_{k,i} & [D_{k,i}-D_{k,0}]F_{k,0}-C_{k,0}\end{bmatrix}
$$

Two observations that make this transparent:

- **$Q_i$ measures the difference between the local and the nominal controller.** $D_{q,i} = D_{k,i} - D_{k,0}$ literally. If $K_i \equiv K_0$ then $Q_i \equiv 0$. The YK parameter is the "correction" the local design applies to the nominal.
- **$A_{q,i}$ is block upper-triangular with two diagonal blocks:**
  - top-left $2\times2$: $\begin{bmatrix}A + B_2D_{k,i}C_2 & B_2C_{k,i}\\ B_{k,i}C_2 & A_{K,i}\end{bmatrix} = A_{cl,i}$ — the **local closed loop**, exponentially stable by (A.2.1);
  - bottom-right: $A_{K,0}+B_{k,0}F_{k,0}$ — the **stabilized nominal controller dynamics**, exponentially stable by LMI (7).

So $Q_i$ is stable *by construction*, with the two LMIs and (A.2.1) supplying exactly the two diagonal blocks needed. No further condition is imposed on the local designs. This is where the decoupling is bought.

### 5.7 The fixed part $J(\rho)$ (eq. 11)

$$
J(\rho) = \left[\begin{array}{cc|cc}
A + B_2F_g & 0 & 0 & B_2\\
-B_{k,0}C_2 & A_{K,0}(\rho,\dot\rho) & B_{k,0} & 0\\ \hline
F_g - D_{k,0}C_2 & C_{k,0} & D_{k,0} & I\\
-C_2 & -F_{k,0} & I & 0
\end{array}\right]
$$

$\tilde K(\rho,\sigma) = \mathcal{F}_l(J(\rho), Q(\rho,\sigma))$ with $Q(\rho,\sigma) = \sum_i \sigma_i(\rho)Q_i(\rho)$. In particular $K_m(\rho) \equiv \mathcal{F}_l(J(\rho),Q_m(\rho))$ when $\sigma_m=1$ and the rest are zero — i.e. the architecture **exactly reproduces** each local controller, it does not approximate it.

$J$ depends only on $G$, $K_0$, $F_g$, $F_{k,0}$. It is **frozen** across all subregions.

### 5.8 Proof of Theorem II.1

**Step 1 — $Q(\rho,\sigma)$ is exponentially stable $\forall\rho, \forall\sigma$.**
Each $A_{q,i}$ is block-triangular with exponentially stable diagonal blocks (§5.6). Lemma A.1 (Appendix, §7 below) then gives a parameter-dependent $X(\rho)$ certifying exponential stability of the full triangular $A_{q,i}$. Since $A_q = \mathrm{diag}(A_{q,i})$ is $\sigma$-independent and block-diagonal with each block stable, $Q(\rho,\sigma) = \sum_i\sigma_iQ_i(\rho)$ is exponentially stable over $\mathcal{P}$ for **every** bounded $\sigma_i\in[0,1]$.

**Step 2 — closed-loop exponential stability.** Apply the state transformation

$$
T = \begin{bmatrix}I&0&0&0\\0&0&0&I\\ I&-I&0&0\\0&0&I&0\end{bmatrix},\qquad
T^{-1} = \begin{bmatrix}I&0&0&0\\ I&0&-I&0\\0&0&0&I\\0&I&0&0\end{bmatrix}
$$

to $CL(\rho,\sigma)$. This yields the **block upper-triangular** form (eq. 13):

$$
\bar A_{cl}(\rho,\sigma) = TA_{cl}T^{-1} =
\begin{bmatrix}
A+B_2F_g & B_2C_q(\rho,\sigma) & \star & \star\\
0 & A_q(\rho,\dot\rho) & B_qC_2 & -B_qF_{k,0}\\
0 & 0 & A+B_2D_{k,0}C_2 & B_2C_{k,0}\\
0 & 0 & B_{k,0}C_2 & A_{K,0}(\rho,\dot\rho)
\end{bmatrix}
$$

**This is the heart of the paper.** The three diagonal blocks are:

| Block | Stability from |
|---|---|
| $A + B_2F_g$ | LMI (6), with $Y_g = X_g^{-1}$ → inequality (14a) |
| $A_q(\rho,\dot\rho)$ | Step 1 → inequality (14b) |
| $A_{cl,0} = \begin{bmatrix}A+B_2D_{k,0}C_2 & B_2C_{k,0}\\ B_{k,0}C_2 & A_{K,0}\end{bmatrix}$ | (A.2.1): $K_0$ stabilizes $G$ → inequality (14c) |

**The switching signal $\sigma$ appears only in the strictly-upper-triangular blocks** ($B_2C_q$ at position (1,2), and the $\star$ entries $-B_2(F_g - (D_{k,0}+D_q)C_2)$ and $B_2(C_{k,0}-D_qF_{k,0})$). Since $\sigma_i\in[0,1]$ is bounded, those off-diagonal terms are bounded, and Lemma A.1 makes bounded off-diagonal terms irrelevant to stability.

> **In one sentence: the YK embedding pushes the switching out of the spectrum of the closed-loop state matrix and into bounded off-diagonal coupling, where it cannot destabilize anything.**

Hence $X_{cl}(\rho) = T^T\,\mathrm{diag}(Y_g, Y_q, Y_{cl,0})\,T$ is a **single, $\sigma$-independent** parameter-dependent Lyapunov certificate valid across all subregions. Compare Lu & Wu [26], where the Lyapunov function *is* switching-dependent and therefore needs junction conditions.

**Performance.** With $V(x_{cl}) = x_{cl}^TX_{cl}(\rho)x_{cl}$ and continuous closed-loop states at every switching instant $t_k$:

$$
V(x_{cl}(t_k)) = V(x_{cl}(t_k^-)) \;\Rightarrow\; V(x_{cl}(t_k)) \le V(x_{cl}(t_k^-)) \tag{15}
$$

— the Lyapunov function does **not jump** at switches, because there is no state reset. (In [26] this inequality must be *enforced* by design constraints; here it is automatic.)

By YK performance recovery, $CL(\rho,\sigma)$ has the performance of $CL_i(\rho)$ inside each $\mathcal{P}_i$. Combining (16) with (15), integrating $\dot V + \frac{1}{\gamma_\infty}z^Tz - \gamma_\infty w^Tw < 0$ from $x_{cl}(0)=0$ and using $V(x_{cl}(T))\ge0$ gives

$$
\|z\|_2 < \gamma_\infty\|w\|_2, \qquad \gamma_\infty = \max_i \{\gamma_{\infty,i}\}
$$

The global level is the **worst local level** — not the level of a single global design. This is the formal statement of "reduced conservatism": you inherit the best of your patchwise designs rather than a compromise across the whole region.

### 5.9 Remark II.3 — what "smooth switching" means here

An important honesty note in the paper: because $\sigma_i$ jumps between $\{0,1\}$, the controller matrices $C_q$ and $D_q$ **do change discontinuously** at a boundary crossing. Smoothness is claimed at the level of **closed-loop states and signals**, not controller matrices. The triangular interconnection guarantees continuous closed-loop states for any bounded measurable switching. That is precisely the property that matters physically (no steering torque spike, no actuator saturation, no mechanical shock).

---

## 6. Contribution II — Partitioned polytopic-based LPV–YK control

This section does something different in kind, not just in technique: it combines **YK gain-scheduling of LTI vertex controllers inside each subset** with **YK switching between subsets**, and the partition is drawn **along the parameter trajectory** to kill over-bounding.

### 6.1 Extra model assumptions

$$
G(\rho):\quad
\begin{cases}
\dot x = A(\rho)x + B_1(\rho)w + B_2u\\
z = C_1(\rho)x + D_{11}(\rho)w + D_{12}u\\
y = C_2x + D_{21}w
\end{cases} \tag{18}
$$

1. strictly proper: $D_{22}(\rho) = 0$;
2. **$B_2$, $C_2$, $D_{12}$, $D_{21}$ parameter-independent** (Angelis [1]).

Assumption 2 is the real restriction versus the grid approach. The paper notes it is not serious in practice: it can always be met by filtering $u$ and $y$ (Apkarian, Gahinet & Becker [3]) or absorbing the parameter dependence via state augmentation.

### 6.2 The partition — where the over-bounding goes away

$\mathcal{P}_0$ is a convex polytope containing all trajectories. The subsets

$$
\mathcal{P}_i := \mathcal{CO}\{w_{i1},\dots,w_{i2^{n_p}}\},\qquad
\rho = \sum_{j=1}^{2^{n_p}}\alpha_{ij}(\rho)w_{ij},\quad \sum_j\alpha_{ij}=1,\ \alpha_{ij}\ge0
$$

are laid out **along the parameter trajectory** and intersect at their **vertices**. For the vehicle, the trajectory is the hyperbola $\rho_2 = 1/\rho_1$; a chain of small polytopes hugging that curve covers vastly less fictitious territory than the single hull $\mathcal{P}_0$. That is the mechanism by which this architecture reduces over-bounding — not by a cleverer LMI, but by **only ever asking the local designs to cover the region the plant actually visits**.

$G(\rho)$ is represented equivalently over $\mathcal{P}_0$ via $\alpha_{0j}$ (eq. 21) or over each $\mathcal{P}_i$ via $\alpha_{ij}$ (eq. 22).

### 6.3 Assumptions (A.3.1) and (A.3.2)

**(A.3.1)** There exists a polytopic LPV controller $K_0(\rho)$ quadratically stabilizing $G(\rho)$ over the **full** $\mathcal{P}_0$ (standard Apkarian–Gahinet–Becker design), realized as $\sum_j \alpha_{0j}(\rho)\begin{bmatrix}A_{K,0j}&B_{k,0j}\\C_{k,0j}&D_{k,0j}\end{bmatrix}$. Again: this may be conservative; it is the net, not the performance.

**(A.3.2) — the boundary-sharing condition, and the cleverest practical idea in the paper.** At each vertex $w_{ij}$, an LTI controller $K_{ij}$ is designed separately to stabilize the LTI plant $G_{ij}$. **At each intersecting boundary a *single, unique* LTI controller is designed** — since $\mathcal{P}_i$ and $\mathcal{P}_{i+1}$ share vertices, one takes

$$
K_{i3} \equiv K_{(i+1)2}
$$

"and not to be designed twice. Even if it is designed twice, it should be done with the same LMI formulation and the same weighting performance."

**Consequence:** at the switching surface, the controller realized from the left is *identical* to the controller realized from the right. The state and control-input energies immediately before and after the switch are unchanged. Switching at a shared face is not merely bounded — it is a **no-op**. This is why the experimental hysteretic chattering between $\mathcal{P}_3$ and $\mathcal{P}_4$ (§10) does no damage.

### 6.4 Lemma III.1 (Bianchi & Peña [10]) — the enabling trick

> For matrices $A_i$ at the vertices of a convex hull $\mathcal{J}$, the following are equivalent:
> 1. $A_i$ is Hurwitz $\forall i$;
> 2. there exist transformation matrices $Z_i$ such that $\bar A(\rho) = \sum_i\alpha_i(\rho)Z_iA_iZ_i^{-1}$ is **quadratically stable** $\forall\rho\in\mathcal{J}$.

Read this carefully, because it is doing a lot of work. Pointwise Hurwitz-ness of vertex matrices does **not** in general imply quadratic stability of their convex combination — this is the classic counterexample territory of LPV theory. Lemma III.1 says: it *does*, **provided you are free to choose the state-space realization at each vertex**.

You cannot use this on a *plant* (its realization is given by physics). But you **can** use it on the **YK parameter**, because $Q$ is a synthetic object whose realization you choose. Similarity transformations leave each vertex transfer function untouched, so $\bar Q_{ij}$ and $Q_{ij}$ are the same system at the vertices — while the *interpolation* of the transformed realizations is a different (and now quadratically stable) LPV system. That is the whole reason the partitioned-polytopic scheme can be certified with a constant Lyapunov matrix.

The paper notes the choice of $Z_{ij}$ (hence $Y_q$) changes the dynamical properties of $\bar Q_\sigma(\rho)$ even though quadratic stability always holds — a design freedom flagged as open in the conclusion.

### 6.5 Theorem III.1 — the synthesis conditions

> Under (A.3.1)–(A.3.2), the switched LPV–YK controller $\tilde K_\sigma(\rho)$ with realizations (31) **quadratically stabilizes** $G(\rho)$ for any $\rho\in\mathcal{P}$, if there exist symmetric positive-definite **constant** matrices $X_g$, $X_{k,0}$, $X_{q,ij}$ and matrices $W_{0j}$, $V_{0j}$ such that $\forall j\in\mathbb{I}[1,2^{n_p}]$:

$$
A_{0j}X_g + X_gA_{0j}^T + B_2W_{0j} + W_{0j}^TB_2^T < 0 \tag{28}
$$
$$
A_{K,0j}X_{k,0} + X_{k,0}A_{K,0j}^T + B_{k,0j}V_{0j} + V_{0j}^TB_{k,0j}^T < 0 \tag{29}
$$
$$
X_{q,ij}A_{q,ij} + A_{q,ij}^TX_{q,ij} < 0 \quad \forall i\in\mathbb{I}_N \tag{30}
$$

with $F_{g,0j} = W_{0j}X_g^{-1}$, $F_{k,0j} = V_{0j}X_{k,0}^{-1}$, $F_{g,ij} = F_{g,0}(w_{ij})$, $F_{k,ij} = F_{k,0}(w_{ij})$, and

$$
A_{q,ij} = \begin{bmatrix}
A_{ij}+B_2D_{k,ij}C_2 & B_2C_{k,ij} & B_2[D_{k,ij}-D_{k,0}(w_{ij})]F_{k,0j}-B_2C_{k,0}(w_{ij})\\
B_{k,ij}C_2 & A_{K,ij} & B_{k,ij}F_{k,0j}\\
0 & 0 & A_{K,0}(w_{ij})+B_{k,0}(w_{ij})F_{k,ij}
\end{bmatrix}
$$

**Structural parallels to the grid case — deliberate and worth noting:**

- (28) ↔ (6): plant state-feedback stabilizability, now **quadratic** (one constant $X_g$ shared across vertices $j$, with per-vertex $W_{0j}$).
- (29) ↔ (7): same for the nominal controller dynamics.
- (30) is **new**: per-vertex, per-subset Lyapunov matrices $X_{q,ij}$ certifying each $A_{q,ij}$ Hurwitz, which then feed Lemma III.1 through $Z_{ij} = X_{q,ij}^{1/2}$.
- $A_{q,ij}$ has exactly the same block-triangular anatomy as $A_{q,i}$ in the grid case: local vertex closed loop on the top-left, stabilized nominal controller dynamics on the bottom-right. LMI (30) is in practice a *verification* rather than a constraint — the paper notes "$A_{q,ij}$ are always Hurwitz by construction (YK concept)".

**Note the index range: (30) ranges over $i\in\mathbb{I}_N$ and all $j$ — this is the only condition whose count grows with the number of subsets, and it is a set of *independent, decoupled* small LMIs, one per vertex. Nothing couples subset $i$ to subset $i'$.** That is the plug-and-play property in the polytopic setting.

### 6.6 The controller realizations (eq. 31)

$$
\tilde A_{K,i}(\rho) = \sum_j\alpha_{ij}(\rho)\begin{bmatrix}
A_{ij}+B_2F_{g,ij}-B_2\bar D_{q,ij}C_2 & -B_2\bar D_{q,ij}F_{k,ij} & B_2\bar C_{q,ij}\\
-B_{k,0}(w_{ij})C_2 & A_{K,0}(w_{ij}) & 0\\
-\bar B_{q,ij}C_2 & -\bar B_{q,ij}F_{k,ij} & \bar A_{q,ij}
\end{bmatrix}
$$

$$
\tilde B_{k,i} = \sum_j\alpha_{ij}\begin{bmatrix}B_2\bar D_{q,ij}\\ B_{k,0}(w_{ij})\\ \bar B_{q,ij}\end{bmatrix},\qquad
\tilde D_{k,i} = \sum_j\alpha_{ij}[D_{k,0}(w_{ij})+\bar D_{q,ij}]
$$

with the **transformed** YK matrices

$$
\bar A_{q,ij} = Z_{ij}A_{q,ij}Z_{ij}^{-1},\quad \bar B_{q,ij}=Z_{ij}B_{q,ij},\quad \bar C_{q,ij}=C_{q,ij}Z_{ij}^{-1},\quad Z_{ij}=X_{q,ij}^{1/2}
$$

and the fixed part $J(\rho) = \sum_j\alpha_{ij}(\rho)\,[\cdot]$ of eq. (33), so that $\tilde K_i(\rho) = \mathcal{F}_l(J(\rho),\bar Q_i(\rho))$ on $\mathcal{P}_i$.

Compare with the grid case: there $\sigma$ selected among $N$ parallel $Q_i$'s; here $\alpha_{ij}$ **convexly interpolates** vertex YK parameters inside the active subset, and $\sigma$ selects which subset's interpolation is live.

### 6.7 Proof of Theorem III.1

**Step 1 — $\bar Q_i(\rho)$ is quadratically stable.** By (30) each $A_{q,ij}$ is Hurwitz (structurally guaranteed by the block-triangular YK construction). Lemma III.1 then supplies $Z_{ij}$ such that $\bar Q_\sigma(\rho)$ is quadratically stable over $\mathcal{P}$. Concretely, with $X_{q,ij} = Z_{ij}^TY_qZ_{ij}$ and the choice $Z_{ij} = X_{q,ij}^{1/2}$ one gets $Y_q = I$:

$$
Y_q\bar A_{q,ij} + \bar A_{q,ij}^TY_q < 0 \quad \forall i,j \qquad\text{with } Y_q = I
$$

A **single constant identity** Lyapunov matrix for all vertices of all subsets. Elegant, and it is what makes the subsequent switching argument free.

**Step 2 — closed-loop quadratic stability.** The same transformation $T$ gives the block-triangular $\bar A_{cl,i}(\rho)$ (eq. 35), with diagonal blocks

| Block | Stability from |
|---|---|
| $A_{ij}+B_2F_{g,ij}$ | (28) with $Y_g=X_g^{-1}$, $W_{0j}(w_{ij})=F_{g,ij}X_g$ → (36) |
| $\bar A_{q,ij}$ | Step 1, $Y_q=I$ → (37) |
| $\mathcal{A}_{ij}=\begin{bmatrix}A_{ij}+B_2D_{k,0}(w_{ij})C_2 & B_2C_{k,0}(w_{ij})\\ B_{k,0}(w_{ij})C_2 & A_{K,0}(w_{ij})\end{bmatrix}$ | (A.3.1): $K_0$ quadratically stabilizes $G$ over $\mathcal{P}_0$ → (38) |

Each of (36)–(38) holds at every vertex and therefore, by convexity, everywhere in the subset. Block-triangularity (Xie & Eisaka [43], Lemma 2) then yields

$$
X_{cl} = T^T\mathrm{diag}(Y_g,Y_q,Y)T,\qquad X_{cl}\tilde A_{cl,i}(\rho)+\tilde A_{cl,i}^T(\rho)X_{cl}<0\quad\forall\rho\in\mathcal{P},\ \forall i\ge1
$$

$X_{cl}$ is **constant and common to every subset $i$**. A switch from $\mathcal{P}_i$ to $\mathcal{P}_{i'}$ therefore requires no junction condition, no dwell time, no hysteresis: the same quadratic Lyapunov function decreases on both sides. Add (A.3.2) — identical controllers on the shared face — and the input–output performance is *conserved exactly* at the transition instant.

---

## 7. Lemma A.1 — the triangular-stability workhorse

Both proofs bottom out here (Rajamani, Nagpal & Khargonekar [32]).

> Let $A(\rho)=\begin{bmatrix}A_{11}(\rho)&A_{12}(\rho)\\0&A_{22}(\rho)\end{bmatrix}$ with $A_{12}$ **bounded**, and $X_1$, $X_2$ bounded positive-definite with
> $$\dot X_1 + X_1A_{11}+A_{11}^TX_1 \le -\alpha_1I,\qquad -\alpha_2I \le \dot X_2+X_2A_{22}+A_{22}^TX_2 < 0$$
> Then there exists $X(\rho)>0$ with $\pi = \dot X + XA + A^TX < 0$.

**Proof sketch.** Choose $X(\rho)=\mathrm{diag}(X_1(\rho), \lambda X_2(\rho))$ and note $\pi$ has $\pi(1,1)\le-\alpha_1I$, $\pi(2,2)=\lambda(\dot X_2+X_2A_{22}+A_{22}^TX_2)<0$ for any $\lambda>0$, and off-diagonal $X_1A_{12}$ bounded by some $\alpha_3$. Schur complement:

$$
\pi \le -\alpha_1I + \lambda^{-1}\alpha_2\alpha_3^2I < 0 \quad\text{for any }\lambda > \alpha_2\alpha_3^2/\alpha_1
$$

**Interpretation, and why it matters twice over.** Stability of a block-triangular system is decided **entirely by its diagonal blocks**; bounded off-diagonal coupling can always be dominated by scaling the Lyapunov weight of the lower block. Used:

- in **Step 1** of both proofs, to conclude each $Q_i$ is stable from its two stable diagonal blocks;
- in **Step 2** of the grid proof, to absorb the $\sigma$-dependent off-diagonal terms of $\bar A_{cl}(\rho,\sigma)$ — which is exactly the "switching cannot destabilize" argument.

The boundedness requirement on $A_{12}$ is why the theorem needs $\sigma_i$ **bounded** and nothing else. No continuity, no dwell time, no rate limit on $\sigma$.

---

## 8. Side-by-side comparison of the two architectures

| | **Grid-based LPV–YK** (Sec. II) | **Partitioned polytopic LPV–YK** (Sec. III) |
|---|---|---|
| Plant class | General LPV; $B_2,C_2,D_{12},D_{21}$ **may vary** with $\rho$ | Strictly proper; $B_2,C_2,D_{12},D_{21}$ **constant** |
| Local controllers | LPV controllers $K_i(\rho)$ over subregions | LTI controllers $K_{ij}$ at subset vertices |
| Certificate | Parameter-dependent $X_{cl}(\rho)$ | **Constant** $X_{cl}$ (quadratic stability) |
| Design conditions | pLMIs (6)–(7) on a grid, rate-bounded | Vertex LMIs (28)–(30) |
| Rate bounds $\dot\rho$ | **Required** ($\underline\nu,\overline\nu$) | **Not required** — a genuine advantage when rates are fast/unknown |
| Offline cost | Higher; scales with grid density | Lighter; scales with vertex count $2^{n_p}$ |
| Conservatism source removed | Single $X(\rho)$ over the wide region | Single $X$ **and** over-bounding (trajectory-aligned partition) |
| Switching mechanism | $\sigma_i$ selects among $N$ parallel $Q_i$ | $\alpha_{ij}$ interpolates within subset; $\sigma$ selects subset |
| Boundary behaviour | Bumpless via triangular structure (Rmk II.3) | **Exactly identical controller** on shared face (A.3.2) |
| Guarantee | Exponential stability + $\gamma_\infty=\max_i\gamma_{\infty,i}$ | Quadratic stability, performance conserved at transitions |
| Plug & play | Add $Q_{N+1}$ to the bank | Add a subset + its vertex $Q$'s |

**Choosing between them.** Grid-based when the plant has parameter-varying input/output matrices, or when you want the lower conservatism of a parameter-dependent certificate and can afford the offline cost. Partitioned polytopic when the parameters are coupled and over-bounding dominates, when scheduling rates are fast or uncertain (no rate bounds needed), or when you want the exact-continuity guarantee at boundaries.

---

## 9. Application: lateral control of an autonomous vehicle

### 9.1 Why this plant

The lateral bicycle model is scheduled by $v_x$ and $1/v_x$ — **coupled, non-affine** parameters. The convex hull over $[5,30]$ m/s contains combinations like $(v_x = 30, 1/v_x = 1/5)$ that are physically impossible. This is precisely the over-bounding pathology. Additionally, prior work [8] showed all three classical LPV approaches became **noise-sensitive above 18 m/s**, traced to the conservatism of a single controller covering up to 30 m/s.

### 9.2 Model

Gridded form (41), states $x = [v_y\ \ \omega]^T$, input $\delta$ (front steering angle):

$$
A_\Sigma(\rho) = \begin{bmatrix}
-\frac{C_r+C_f}{m\rho} & -\frac{C_fl_f-C_rl_r}{m\rho}-\rho\\[4pt]
-\frac{C_fl_f-l_rC_r}{I\rho} & -\frac{C_fl_f^2+l_r^2C_r}{I\rho}
\end{bmatrix},\qquad
B_\Sigma = \begin{bmatrix}\frac{C_f}{m}\\ \frac{C_fl_f}{I}\end{bmatrix},\qquad
C_\Sigma=[0\ \ 1]
$$

Affine/convex form (42), $\rho_1 = v_x$, $\rho_2 = 1/v_x$:

$$
A_\Sigma(\rho) = \begin{bmatrix}
-\frac{C_r+C_f}{m}\rho_2 & -\frac{C_fl_f-C_rl_r}{m}\rho_2 - \rho_1\\[4pt]
-\frac{C_fl_f-l_rC_r}{I}\rho_2 & -\frac{C_fl_f^2+l_r^2C_r}{I}\rho_2
\end{bmatrix}
$$

(41) is used for the grid designs, (42) for the polytopic ones. The model is augmented with an identified steering actuator $\Sigma_{act}$.

### 9.3 Design setup

- $\rho \in \mathcal{P} = [5,30]$ m/s, $|\dot\rho| \le 5$ m/s².
- Five subsets: $[5,10]\cup[10,15]\cup[15,20]\cup[20,25]\cup[25,30]$.
- Weights — **deliberately different for nominal vs local**:
  - nominal: $W_{e,0} = \frac{q+2}{2q+0.002}$, $W_{u,0} = \frac{q+5}{0.01q+5}$
  - local: $W_e = \frac{q+2}{2q+0.002}$, $W_u = \frac{q+10}{0.01q+10}$

  The local designs get a more permissive $W_u$ (higher corner frequency) — i.e. more control authority — because they only need to hold a narrow speed band. The nominal is detuned, as befits a global safety net.

**Grid design steps:** design $K_0$ over $\mathcal{P}$ with $(W_{e,0},W_{u,0})$ → $\gamma_{\infty,0}=1.55$; design $K_i$ over each $\mathcal{P}_i$ with $(W_e,W_u)$; solve (6)–(7) with the **affine** Lyapunov basis

$$
X_g(\rho) = X_g^0 + X_g^1\rho,\qquad X_{k,0}(\rho)=X_{k,0}^0+X_{k,0}^1\rho
$$

→ $F_g(\rho)$, $F_{k,0}(\rho)$; build $J(\rho)$ from (11) and $Q_i(\rho)$ from (10); set $Q(\rho,\sigma)=\sum_{i=1}^5\sigma_i(\rho)Q_i(\rho)$ with $\sigma_i\in\{0,1\}$ toggled at boundaries.

**Table I — local $\gamma$-performances (grid):**

| $\gamma_{\infty,1}$ | $\gamma_{\infty,2}$ | $\gamma_{\infty,3}$ | $\gamma_{\infty,4}$ | $\gamma_{\infty,5}$ |
|---|---|---|---|---|
| 1.1543 | 1.1157 | 1.2150 | 1.2112 | 1.2175 |

A single grid-based LPV controller over the full region with the same weights achieves $\gamma_\infty = 1.42$. So $\max_i\gamma_{\infty,i} = 1.2175$ versus $1.42$ — **a ~14% improvement in guaranteed $\mathcal{L}_2$ level**, and this is the quantitative payoff of Theorem II.1's $\gamma_\infty = \max_i\gamma_{\infty,i}$ result.

**Polytopic design steps:** $\mathcal{P}_0$ is the **triangle** $\mathcal{CO}\{(\underline\rho_1,\underline\rho_2),(\underline\rho_1,\overline\rho_2),(\overline\rho_1,\underline\rho_2)\}$ — a triangle rather than a rectangle already trims the unreachable corner. Five triangular subsets along the $\rho_2=1/\rho_1$ trajectory. LTI $K_{ij}$ at all vertices with $(W_e,W_u)$; polytopic $K_0$ with $(W_{e,0},W_{u,0})$ → $\gamma_{\infty,0}=1.6$. $F_g(\rho)$ and $F_{k,ij}$ via LMI state-feedback (pole-placement constraints or LQR). $J(\rho)$, $\bar Q(\rho,\sigma)$ from (33), (32).

**Table II — LTI $\gamma$-performances:**

| | $w_{i1}$ | $w_{i2}$ | $w_{i3}$ |
|---|---|---|---|
| $K_{1j}$ | 0.99 | 1.12 | 1.00 |
| $K_{2j}$ | 0.97 | 1.00 | 0.98 |
| $K_{3j}$ | 0.97 | 0.98 | 0.98 |
| $K_{4j}$ | 0.97 | 0.98 | 0.98 |
| $K_{5j}$ | 0.98 | 0.98 | 0.99 |

All close to ~0.98 versus $\gamma_{\infty,0}=1.6$ for the global polytopic design — a striking illustration of how much the over-bounded hull was costing.

---

## 10. Experimental results

**Platform.** Robotized electric Renault ZOE, computer-controlled steering and pedals, GPS + IMU, dSPACE MicroAutoBox, private test track at Satory (France) — including poor road surface and inclinations, so robustness is genuinely exercised. A navigation restart at the end of the highway section produces the excluded "NAV" band in the plots.

**Benchmarks.** (1) the Lu & Wu [26] switched-LPV controller; (2) a gain-scheduled LPV–YK controller built by scheduling LTI controllers at the vertices of the **full** $\mathcal{P}_0$ (in the style of [13]) — i.e. YK but *without* partitioning, which isolates the contribution of the partition itself.

**Grid-based LPV–YK vs switched LPV [26]** (Fig. 7). Both minimize lateral error comparably (Fig. 7a) — as expected, since both use similar $H_\infty$ design principles. The difference is at the switching instants: zoomed views of yaw rate and steering angle (Figs. 7b, 7c) show the switching effect is **negligible** for the LPV–YK controller and clearly visible for the switched-LPV controller. Control input at high speed stays bounded and consistent with the designed bandwidth — the high-speed noise sensitivity reported in [8] is resolved.

**Partitioned polytopic LPV–YK vs gain-scheduled LPV–YK** (Fig. 8). Both achieve lateral error minimization (Fig. 8a). The gain-scheduled version — same YK machinery, no partition — produces **markedly higher steering oscillations** (Figs. 8b, 8c). This is the clean experimental isolation of the partitioning benefit: YK alone is not enough; the trajectory-aligned partition is what removes the over-bounding.

**The hysteresis episode.** At $t\in[28.3, 28.5]$ s the scheduling signal chatters back and forth between $\mathcal{P}_3$ and $\mathcal{P}_4$. This would be a worst case for a dwell-time-based scheme. Here it has **no effect** on steering input or yaw rate, because at that shared boundary

$$
\mathcal{F}_l(G, J_{33}, Q_{33}) \equiv \mathcal{F}_l(G, J_{42}, Q_{42})
$$

the two controllers are dynamically identical by (A.3.2). This is the single most convincing experimental validation in the paper: an adversarial switching pattern, deliberately chattering, with a null effect.

**The honest caveat (Sec. IV-D).** When the vehicle operates **well inside** a subset with no transition, the proposed controllers perform *comparably* to the benchmarks — all rely on similar LPV/$H_\infty$ principles. The benefit appears **near boundaries**. The paper states this plainly rather than overclaiming.

**Claimed advantages, restated with their basis:**

| Advantage | Basis |
|---|---|
| Stability by construction | YK preserves each pre-designed controller's stability; no re-verification of global/multiple Lyapunov functions |
| Smooth switching | Triangular interconnection + (A.3.2) identical boundary controllers; validated experimentally |
| Reduced conservatism | Lyapunov function required per partition, not over the whole region as in [10] |
| Plug & play | Add $Q_k(\rho)$; stability and continuity follow without altering existing $Q_j$'s — vs re-certifying dwell-time logic or recomputing global Lyapunov functions in [10], [16], [26] |

**On the assumptions' practicality (Sec. IV-D).** Most plants are strictly proper; parameter-dependent $B_2,C_2,D_{12},D_{21}$ can be absorbed by state-augmentation preprocessing; operating sets are compact by design; scheduling rates are usually slower than the control bandwidth, and when they are fast or uncertain, the polytopic formulation **avoids rate bounds entirely**. A global nominal $K_0(\rho)$ is routinely obtainable by robust LPV/$H_\infty$ synthesis, possibly conservative — with local performance recovered through the YK parameters.

---

## 11. Mapping theory to the code in this repository

Every equation in the theory above has a named function in `src/`. The two `examples/` scripts run the pipelines end to end.

### 11.1 Equation-to-function map

| Paper | Function | Location |
|---|---|---|
| Plant (41), gridded non-affine form | `bicycleModelGrid` | `src/model/` |
| Plant (42), affine form $\rho_1=v_x$, $\rho_2=1/v_x$ | `bicycleModelPolytopic` | `src/model/` |
| $\Sigma_{act}$, 2nd-order + Padé(2) delay | `steeringActuator` | `src/model/` |
| Weights $W_e$, $W_u$ (44)–(45) | `weightPerformance`, `weightControl` | `src/weights/` |
| Partition (43) | `partitionInterval` | `src/utils/` |
| Triangular $\mathcal{P}_0$ (Fig. 6) | `polytopeVertices` | `src/utils/` |
| $\sum_j \alpha_{ij}(\rho)(\cdot)$ (20) | `polytopicInterpolation` | `src/utils/` |
| $K_0(\rho)$, $K_i(\rho)$ grid designs | `designGridLpvController` | `src/synthesis/` |
| $K_0(\rho)$ polytopic, (A.3.1) | `designPolytopicLpvController` | `src/synthesis/` |
| $K_{ij}$ LTI vertex designs, (A.3.2) | `designLtiController` | `src/synthesis/` |
| LMIs (6)–(7) → $F_g$, $F_{k,0}$ | `lmiStateFeedbackGrid` | `src/lmi/` |
| LMIs (28)–(29) → $F_{g,0j}$, $F_{k,0j}$ | `lmiHinfStateFeedbackPolytope` | `src/lmi/` |
| Region-constrained state feedback | `lmiStateFeedbackRegion` | `src/lmi/` |
| Polytopic LPV H∞ synthesis | `lmiHinfPolytope` | `src/lmi/` |
| $J(\rho)$ (11), (33) via coprime factors | `coprimeFactorizationLpv`, `coprimeFactorization` | `src/youla/` |
| $Q_i(\rho)$ (10), (32) | `youlaParameter` | `src/youla/` |
| $\tilde K = (U+MQ)(V+NQ)^{-1}$ (8), (31) | `ykController` | `src/youla/` |
| Fixed $n_q$ for $A_q=\mathrm{diag}(A_{q,i})$ (9) | `equalizeOrder` | `src/utils/` |
| Factorization identity checks | `verifyStepEquality` | `src/utils/` |

### 11.2 Notes on the implementation choices

- **Coprime-factor realization.** The controller is realized as $(U_0+M_0Q)(V_0+N_0Q)^{-1}$ rather than the explicit LFT state-space of (8)/(31). These are equivalent, but the factor form makes the architectural claim concrete: $M_0, N_0, U_0, V_0^{-1}$ are *fixed* blocks and only $Q$ is switched. In a Simulink/dSPACE realization that means the switching logic touches one subsystem.

- **`youlaParameter` takes the nominal controller evaluated at the local operating point.** In the polytopic pipeline this is $K_0(w_{ij})$, obtained by interpolating the nominal vertex controllers at $w_{ij}$ — matching the $A_{K,0}(w_{ij})$, $C_{k,0}(w_{ij})$, $D_{k,0}(w_{ij})$ terms throughout (31)–(32).

- **`equalizeOrder` pads with fast decoupled modes.** Locally designed $Q_i$ generally have different McMillan degree, but the block-diagonal bank needs a constant state dimension. Padding adds modes at $-10^4$ with zero rows in $B$ and zero columns in $C$: uncontrollable, unobservable, and numerically negligible at a 100 Hz sample rate.

- **Shared-face controllers are enforced, not hoped for.** `example02` keys vertex designs on their $(\rho_1,\rho_2)$ coordinates and designs each unique vertex exactly once, so $K_{i,3}$ and $K_{i+1,2}$ are literally the same object. It then asserts the identity. This makes (A.3.2) a structural property of the code rather than a consequence of `hinfsyn` being deterministic.

- **Discretization** is `c2d(..., 0.01, 'Tustin')` throughout — 100 Hz, matching the MicroAutoBox target of the original experiments.

- **Loop sign convention: positive feedback.** The `sysic` interconnection in every synthesis function defines the measurement as `-r + G + n`, i.e. the controller is fed $y - r$ and closes the loop as $u = K\,y$ when $r = 0$. The coprime factorization uses the same convention — its well-posedness terms are $(I - D_{k}D_{g})^{-1}$, not $(I + D_kD_g)^{-1}$, and the local closed loop inside $A_Q$ is $\begin{bmatrix}A + BD_1C & BC_1\\ B_1C & A_1\end{bmatrix}$. Consequently any closed loop formed outside the LFT must use `feedback(G*K, 1, +1)`; the default `feedback(G*K, 1)` closes the *negative* loop, which is a different (generally unstable) system. The examples and unit tests do this; keep it in mind when adding analysis code.

- **The pole-placement fallback in `coprimeFactorization` is only admissible for stable pairs.** When `Fg` / `Fk0` are passed as `[]`, the gains are chosen to keep the *existing* poles of $G$ and $K_0$. That yields factors in $RH_\infty$ only if $G$ and $K_0$ are already stable, which holds for the vehicle model at the speeds considered but not in general (H-infinity controllers for unstable plants are frequently unstable themselves). The LPV version deliberately provides no fallback; the unit tests use an unstable surrogate precisely so that this requirement is exercised.

### 11.3 Differences from the original research code

This repository is a cleaned reimplementation. The behaviour-relevant differences:

| Change | Rationale |
|---|---|
| Four near-identical `Design_*GridLPV_Controller` functions merged into `designGridLpvController` | They differed only in the weights, which are now arguments |
| Coprime factors returned as a **struct** rather than 8 positional outputs | Call sites become readable and order-independent |
| `hinfsyn` output capture corrected | The original bound its 2nd output to a variable named `GAM`, but `hinfsyn` returns `[K, CL, GAM, INFO]` — so that variable held the *closed loop*, not gamma |
| Gravity constant `9.98` → `9.81` | Typo in the original; unused by the dynamics, so no numerical effect |
| `Weight_e2.m` / `Weight_u2.m` filename-vs-function mismatch removed | Both declared `function ... = Weight_e(...)`; MATLAB resolved by filename, but it is a trap |
| Diagnostic `Create_Glist` / `check_polyt_list` calls removed from the design path | They printed rank checks and returned unused outputs |
| Parameter-varying-weight design variant dropped | `Design_LPVWeights_GRIDController` is an alternate design not part of the published method |
| Benchmarks, Simulink models, experimental data, and RL extensions excluded | Out of scope for a design-pipeline repository |

---

## 12. Critical reading: strengths, limits, editorial notes

### 12.1 What is genuinely strong

1. **The triangularization argument is clean and decisive.** Pushing $\sigma$ into bounded off-diagonal blocks, then invoking Lemma A.1, is a short proof of a strong statement (arbitrary bounded switching). It is not a refinement of dwell-time analysis; it sidesteps the question.
2. **The decoupling is real, not cosmetic.** LMIs (6)–(7) genuinely do not involve $N$, $\sigma$, or $K_i$. The plug-and-play claim follows immediately rather than needing a separate argument.
3. **Assumption (A.3.2) is a small idea with a large effect.** Designing one controller at a shared vertex instead of two converts "bounded transient at the switch" into "no transient at the switch". The hysteresis episode at $t\in[28.3,28.5]$ s is the ideal experiment for it.
4. **Using Lemma III.1 on $Q$ rather than on the plant** is the right insight: realization freedom exists precisely where the object is synthetic.
5. **Honest reporting.** Sec. IV-D states that away from boundaries performance is merely comparable to benchmarks. Remark II.3 explicitly disclaims continuity of controller matrices. Both are the kind of qualification that strengthens credibility.

### 12.2 Limits worth keeping in view

1. **Controller order.** $n_x + n_{k,0} + \sum_{i=1}^N n_{q,i}$, and all $N$ YK parameters run continuously. The code's `-10000*eye(...)` padding, which equalizes $n_q$ across subsets, makes this concrete: order is set by the *largest* $Q_i$, times $N$. Acknowledged in the paper; a real constraint on embedded targets.
2. **Performance level at intermediate $\sigma$.** The bound $\gamma_\infty = \max_i\gamma_{\infty,i}$ rests on YK performance recovery, which is exact when $\sigma$ sits at a vertex of $\{0,1\}^N$ (then $\tilde K \equiv K_l$). For $\sigma$ in the interior — which the *stability* theorem explicitly permits — the closed loop is not equal to any single $CL_i$, and the performance claim is not established with the same rigour. Stability is unconditional; the $\gamma$ bound is cleanest for selector-type $\sigma$. The vehicle implementation uses $\sigma_i\in\{0,1\}$, so this gap is not exercised experimentally.
3. **A small tension between Remark II.3 and Sec. IV-D.** Remark II.3 says smoothness means bumpless *closed-loop states*, and that $C_q$, $D_q$ do jump. The Discussion bullet says "implementing $Q(\rho)$ with continuous partition regions yields continuous $K(\rho)$". Remark II.3 is the accurate statement for the grid architecture with $\sigma_i\in\{0,1\}$; the Discussion sentence reads as describing the polytopic architecture under (A.3.2). Worth disambiguating if revising.
4. **$K_0$ must exist over the whole region.** Both theorems need a globally stabilizing nominal controller. If the region is so wide that *no* single LPV controller stabilizes it, neither architecture applies as stated. The paper's answer (Sec. IV-D) is that a possibly-very-conservative $K_0$ from robust LPV/$H_\infty$ synthesis suffices, since performance is recovered locally — reasonable, but it is a genuine precondition.
5. **Partition design is left to the engineer.** How many subsets, and where to place boundaries, is chosen by inspection (five equal 5 m/s bands). Work like Zhao & Nagamune [44] optimizes switching surfaces; combining that with this framework is an obvious opening.
6. **Choice of $Z_{ij}$ / $F_g$, $F_{k,0}$ is unexploited freedom.** The paper's own conclusion names this as future work: these gains shape the $Q$ dynamics and hence transient behaviour, while stability holds for any admissible choice. The code's two parallel paths for $F_g$ (`lpvsfsyn` vs `lmiStateFeedbackGrid_2016`, with one commented out) suggests this was felt in practice.

### 12.3 Editorial / typographic observations in the manuscript

Useful if this version goes through another revision round:

- Sec. I-C, advantage 2: "without requiring **s** stability proof" — stray character; should read "a stability proof".
- Eq. (7) writes $A_{K,0}(\rho)$ while the controller definition (2) uses $A_{K,i}(\rho,\dot\rho)$; the $\dot\rho$ dependence should be stated consistently (it reappears in (13)).
- After (9): "Notice that $X_g^{-1}(\rho)$ and $X_{k,0}^{-1}(\rho)$ are written smooth with respect to the parameter changes such as $X_g(\rho) = \rho.X_g$" — this sentence is garbled and its intent (that the inverse must remain a smooth, differentiable function of $\rho$, e.g. via an affine basis) is not clear as written.
- Theorem III.1: "having the **stat-space** realizations" → "state-space".
- Sec. III-C: $V_{0j}\in\mathbb{R}^{m\times n_{k,ij}}$ — given $F_{k,0j}=V_{0j}X_{k,0}^{-1}$ with $X_{k,0}\in\mathbb{R}^{n_{k,0}\times n_{k,0}}$, the second dimension should be $n_{k,0}$.
- (A.3.2) is stated for $i\in\mathcal{I}_{N-1}$ while Fig. 3 and the design section use $i\in\mathbb{I}[1,N]$; the $N-1$ presumably refers to the count of *shared boundaries*, but as written it reads as the vertex-design range.
- Eq. (35) contains both $\bar D_{q,ij}$ and $\tilde D_{q,ij}$ in adjacent entries — one appears to be a typo for the other.
- Sec. IV-C and IV-B2 contain unresolved cross-references: "shown in Fig. **??**" appears six times (vehicle photo, test track, speed profile, operating-subset plot). Table II is captioned $i\in[1,5], j\in[1,3]$ while the surrounding text says $\forall i\in\mathbb{I}[1,3], \forall j\in\mathbb{I}[1,5]$ — the indices are transposed in one of the two.
- Sec. IV-B2 says "the LTI state-feedback controllers $F_{k,ij}$, $\forall i\in\mathbb{I}[1,4]$" — should presumably be $\mathbb{I}[1,5]$ to cover all five subsets.
- Remark II.3: "$\sigma_i \in 0,1$" should be $\sigma_i\in\{0,1\}$.
- Fig. 7 and Fig. 8 subfigures (b) and (c) have captions but the rendered plots for (b)/(c) are missing axis labels in the extracted text; confirm they survive typesetting.

---

## Quick reference: the logical skeleton

**Grid-based (Theorem II.1)**

```
(A.2.1) each K_i stabilizes G on P_i        ─┐
LMI (7): A_{K,0} stabilizable by F_{k,0}    ─┼─→ A_{q,i} block-triangular, stable diagonal
                                             │    ──(Lemma A.1)──→ Q_i stable
                                             │    ──→ Q(ρ,σ)=Σσ_i Q_i stable ∀σ bounded
LMI (6): A stabilizable by F_g              ─┤
(A.2.1) K_0 stabilizes G on P               ─┘
                    │
                    ├─ transform T ─→ Ā_cl block-triangular
                    │   diag: (A+B₂F_g), A_q, A_cl,0   ← all stable
                    │   σ lives ONLY in bounded off-diagonal blocks
                    │
                    └─→ single σ-independent X_cl(ρ)
                        ⇒ exponential stability ∀ bounded σ
                        ⇒ ‖z‖₂ < γ_∞‖w‖₂,  γ_∞ = max_i γ_∞,i
```

**Partitioned polytopic (Theorem III.1)**

```
(A.3.2) LTI K_ij at each vertex (shared at boundaries)  ─┐
LMI (29): A_{K,0j} stabilizable by F_{k,0j}             ─┼─→ A_{q,ij} Hurwitz  [LMI (30)]
                                                         │   ──(Lemma III.1, Z_ij=X_{q,ij}^{1/2})──→
                                                         │      Q̄_σ(ρ) quadratically stable, Y_q = I
LMI (28): A_{0j} stabilizable by F_{g,0j}               ─┤
(A.3.1) K_0 quadratically stabilizes G on P_0           ─┘
                    │
                    ├─ transform T ─→ Ā_cl,i block-triangular
                    │   diag: (A_ij+B₂F_g,ij), Ā_q,ij, 𝒜_ij   ← (36),(37),(38)
                    │
                    └─→ constant X_cl common to ALL subsets
                        ⇒ quadratic stability ∀ρ∈P, ∀σ
                        ⇒ (A.3.2) ⇒ identical controller at shared face
                                   ⇒ performance conserved exactly at transitions
```
