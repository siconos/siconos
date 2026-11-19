# Projected Moreau–Jean time-stepping schemes

This note describes the two projected Moreau–Jean schemes of the kernel, their
variants and options, how they map onto the code, and how they were validated.

- `TimeSteppingDirectProjection` + `MoreauJeanDirectProjectionOSI`: **Algorithm 1**
  (direct projection);
- `TimeSteppingCombinedProjection` + `MoreauJeanCombinedProjectionOSI`: **Algorithm 2**
  (combined projection);
- both use `MLCPProjectOnConstraints` for the position projection and the helpers of
  `ProjectionOnConstraintsTools.hpp` for the position update.

Reference: V. Acary, *Projected event-capturing time-stepping schemes for nonsmooth
mechanical systems with unilateral contact and Coulomb's friction*, Computer Methods in
Applied Mechanics and Engineering 256 (2013) 224–250 (below: **CMAME 2013**).

---

## 1. Problem and notation

Lagrangian system with unilateral constraints and Newton impact law:

```
M(q) dv = F(t, q, v) dt + G(q)ᵀ di,     q̇ = v,
0 ≤ g(q) ⊥ λ ≥ 0,                         (position level)
G(q) v⁺ = −e G(q) v⁻   if g(q) = 0         (impact law)
```

| symbol | meaning | in Siconos |
|---|---|---|
| `q, v` | positions, velocities | `ds->q()`, `ds->velocity()` |
| `q_k, v_k` | values at the beginning of the step | DS memories `qMemory()`, `velocityMemory()` |
| `g(q)` | gaps (constraints `g ≥ 0`) | interaction output `y(0)` |
| `G(q) = ∇g(q)ᵀ` | Jacobian of the constraints | `LagrangianR::H()` (computed by `computeJach`) |
| `ġ = G v` | relative velocity | output `y(1)`; at the beginning of the step: memory `y_k(1)` |
| `W` | iteration matrix `M + hθC + h²θ²K` | `MoreauJeanOSI` (`getOSIMatrix`) |
| `θ` | θ-method parameter | `MoreauJeanOSI(theta)` |
| `γ` | forecast parameter of the activation rule | `MoreauJeanOSI::setGamma()`, default 1/2 |
| `e` | restitution coefficient | `NewtonImpactNSL` |
| `h` | time step | `TimeDiscretisation` |
| `λ` (level 1), `μ` (level 0) | velocity impulse, projection multiplier | `lambda(1)`, `lambda(0)` |
| `p(1) = Gᵀλ`, `p(0) = Gᵀμ` | inputs of level 1 and 0 | `ds->p(1)`, `ds->p(0)` |

### Why project?

Moreau–Jean enforces the constraints at the velocity level only:
`0 ≤ G v_{k+1} + e G v_k ⊥ λ ≥ 0` on the contacts considered active. Between an
activation and the impact, the position update
`q_{k+1} = q_k + h(θ v_{k+1} + (1−θ) v_k)` can go below `g = 0` (penetration of
order `h·|v|`), and with persistent contacts the penetration drifts. The projected
schemes add a projection of the positions on the admissible set at the end of the
step. The velocities are **not** modified by the projection.

---

## 2. Algorithm 1: direct projection

Classes: `TimeSteppingDirectProjection`, `MoreauJeanDirectProjectionOSI`.

One step from `(q_k, v_k)`:

1. **Activation** (`MoreauJeanDirectProjectionOSI::addInteractionInIndexSet`, indexSet 1;
   rule `FIRSTWAY_ACTIVATION` compiled in `MoreauJeanDirectProjectionOSI.cpp`). With the
   linear forecast of the gap computed from the outputs at the beginning of the step,

   ```
   g̃_α = g_α(q_k) + γ h ġ_α(q_k, v_k)       (y(0) + γ h y(1))
   ```

   - add α if `g̃_α ≤ activateYPosThreshold()` (default `1e-7`);
   - remove α if `g̃_α > deactivateYPosThreshold()` (1e-7) **and**
     `ġ_α ≥ deactivateYVelThreshold()` (0).

   `γ = 1/2` unless `setGamma()` (or the `(theta, gamma)` constructor) is used.
2. **Velocity step** (`TimeStepping::newtonSolve`, OSNS `SICONOS_OSNSP_TS_VELOCITY`, input/output
   level 1, on indexSet 1):

   ```
   v_free = v_k + h W⁻¹ (θ F_{k+1} + (1−θ) F_k − ...)
   v_{k+1} = v_free + W⁻¹ G_αᵀ λ,   0 ≤ G_α (v_{k+1} + e v_k) ⊥ λ ≥ 0,   G = G(q_k)
   a = q_k + h (θ v_{k+1} + (1−θ) v_k)                       (DS positions after the step)
   ```
3. **Projection** (OSNS `SICONOS_OSNSP_TS_POS` = `MLCPProjectOnConstraints`, input/output
   level 0) on the interactions of `indexSet(_indexSetLevelForProjection)`:
   - level 0 (constructor argument `level = 0`): all the constraints, the setting of
     Algorithm 1;
   - level 1 (default): only the contacts active in the velocity step.

   Loop while a constraint is violated
   (`max(−g_α) > _constraintTolUnilateral` for unilateral constraints,
   `|g_α| > _constraintTol` for equality constraints, both 1e-8) and at most
   `_projectionMaxIteration` (10) times: solve the linearized projection
   (a complementarity problem for unilateral constraints) and update the
   positions with the rule of §4 (`ProjectionUpdate`).
4. The velocities `v_{k+1}` are kept.

Options: `setDoProj(false)` (no projection), `setDoOnlyProj(true)` (projection without
velocity step), `setProjectionMaxIteration`, `setConstraintTol`,
`setConstraintTolUnilateral`, `setProjectionUpdate`, `setProjectionTolq`.

---

## 3. Algorithm 2: combined projection

Classes: `TimeSteppingCombinedProjection`, `MoreauJeanCombinedProjectionOSI`.

The active set is not guessed by a forecast. It is built inside the step: start with
the empty set, do the step, add the constraints violated at the end of the step, and
restart the step with the enlarged set. The velocity problem and the position
projection use **the same active set**, and the projection is an **equality**
projection (the active contacts are held closed: `_doCombinedProjOnEquality = true`).

One step from `(q_k, v_k)`, with `A_0 = ∅` (indexSet 1 and 2 are cleared at the beginning
of `advanceToEvent()`):

```
for i = 0, 1, ...  (index set iterations, _nbIndexSetsIteration ≤ _kIndexSetMax = 50)
  1. velocity step on A_i:
       0 ≤ G_α(q_k) (v_{k+1} + e v_k) ⊥ λ_α ≥ 0,  α ∈ A_i
       a = q_k + h (θ v_{k+1} + (1−θ) v_k)
  2. projection on A_i (equality):
       while |g_α(q_j)| ≥ _constraintTolUnilateral for some α ∈ A_i  (≤ _projectionMaxIteration)
         solve  G_j G_jᵀ μ = −(right-hand side of §4),   update q (§4)
  3. A_{i+1} = A_i ∪ { α : g_α(q_{k+1}) < −_constraintTolUnilateral }
     if A_{i+1} = A_i : stop (index sets stable)
```

Call sequence in `TimeSteppingCombinedProjection::advanceToEvent()`:

```
initialize(); clear indexSet 1 and 2
while (!_isIndexSetsStable):
    _nbIndexSetsIteration++
    TimeStepping::newtonSolve()            # step 1
      └─ initializeNewtonSolve()
           └─ updateAndSwapAllOutput()     # overridden, see below
              updateIndexSets(); initializeOneStepNSProblem(); ...
    updateIndexSet(2)                      # indexSet 2 = indexSet 1
    computeCriteria()                      # y(0) = g(q), H = G(q) on indexSet 2
    [from_origin: saveProjectionOrigin()]
    while (running projection):            # step 2
        [from_origin: shiftProjectionRhs()]
        computeOneStepNSProblem(SICONOS_OSNSP_TS_POS); updateInput(0)
        updatePositionsAfterProjection(); updateWorldFromDS(); computeCriteria()
    update forces of the DS
    updateIndexSet(1)                      # step 3: add violated constraints
```

Points that matter for the correctness of the algorithm:

- **Each index set iteration restarts the step.** `newtonSolve()` begins with
  `initializeNewtonSolve()` → `updateAndSwapAllOutput()`, which computes the outputs at
  the current state and saves them in memory (`y_k`). That is right for the first
  iteration only. For `i ≥ 1`, `TimeSteppingCombinedProjection::updateAndSwapAllOutput()`
  resets the DS states to `(q_k, v_k)` (memories) and recomputes the outputs and
  Jacobians there, **without** swapping the memories. Otherwise (behaviour before the
  fix):
  - `y_k(1)` held the end-of-step relative velocity of the previous iteration, so the
    Newton law used `e·G v_free` instead of `e·G v_k` (bouncing ball:
    `v⁺ = 4.0172` instead of `0.9 × 4.4145 = 3.97305`);
  - the velocity problem was built with `G` at the end of the previous iteration instead
    of `G(q_k)` (rocking block: θ̇ error 2e-2 at the second impact).
- **Only violated constraints are activated** (`MoreauJeanCombinedProjectionOSI::
  addInteractionInIndexSet(inter, 1)`: `g_α < −ε`, `ε = constraintTolUnilateral()`).
  With the former rule `g_α ≤ 0`, a contact closed by the projection (`g = 0` exactly)
  was activated at the next step, and the equality projection then held it closed:
  the bouncing ball stayed glued to the ground (`q = 0` with `v > 0`, `μ < 0`).
- indexSet 2 is the whole indexSet 1 (`addInteractionInIndexSet(inter, 2)` returns
  `true`).

Only indexSet level 2 for the projection and `MoreauJeanCombinedProjectionOSI` are
supported.

---

## 4. The position projection and its variants

### 4.1 The projection problem

Let `a` be the positions after the velocity step. The projection is

```
min_q ½ ‖q − a‖²   s.t.   g_α(q) = 0, α ∈ A     (Algorithm 2, equality on the active set)
                    or    g(q) ≥ 0               (Algorithm 1)
```

with first-order (KKT) conditions

```
q = a + G(q)ᵀ μ,    g_α(q) = 0                     (Algorithm 2)
q = a + G(q)ᵀ μ,    0 ≤ g(q) ⊥ μ ≥ 0               (Algorithm 1)
```

The constraints are nonlinear, so it is solved by a sequence of linearized projections.
Each iteration solves one MLCP (`MLCPProjectOnConstraints`), with matrix `G_j G_jᵀ`,
`G_j = G(q_j)`, `q_0 = a`. The MLCP right-hand side is the output `y(0)`
(`MLCPProjectOnConstraints::computeqBlock`). The solution `μ` is stored in `lambda(0)`, and
`updateInput(0)` computes `p(0) = G_jᵀ μ`.

### 4.2 Variants (`enum class ProjectionUpdate`)

| variant | right-hand side `y(0)` | update | converges to |
|---|---|---|---|
| single linearization (`_projectionMaxIteration = 1`) | `g(a)` | `q = a + G(a)ᵀ μ` | `g(q) = O(κ‖μ‖²)`: not exactly on the manifold |
| **`incremental`** (default) | `g(q_j)` | `q_{j+1} = q_j + G_jᵀ λ_j` | a point of `g = 0` (Gauss–Newton on `g`), **not** the projection |
| **`from_origin`** | `g(q_j) + G_j (a − q_j)` | `q_{j+1} = a + G_jᵀ μ_j` | the fixed point `q = a + G(q)ᵀ μ`, `g = 0`: the projection |

(`κ` = curvature of `g`.)

**incremental.** Each iteration projects the current point `q_j` on the linearized
manifold at `q_j`. The limit `q*` satisfies `g(q*) = 0`, but
`q* − a = Σ_j G(q_j)ᵀ λ_j` is a combination of normals at different points. It is not
normal to the manifold at `q*`, so `q*` is not the closest point to `a`. The distance
between the two limits is `O(κ‖μ‖²)`: second order, but it changes the trajectory,
and on chaotic systems (the rocking block) it is amplified.

**from_origin.** Each iteration linearizes the constraint at `q_j` but projects
**from `a`**:
`g(q_j) + G_j (q − q_j) = 0` with `q = a + G_jᵀ μ` gives
`G_j G_jᵀ μ = −[g(q_j) + G_j (a − q_j)]`.
At a fixed point, `q = a + G(q)ᵀ μ` and `g(q) = 0`: these are exactly the KKT
conditions, so `q` is the projection of `a`. This is the projection of the Scilab code
Hots (`ProjectedMoreau.sci`, `CombinedProjectedMoreau.sci`), which reproduces the
figures of CMAME 2013.

**Linear constraints** (`G` constant): both variants coincide, and one iteration is
enough (ball, column of beads: the results are identical to round-off).

**Cost.** Same matrix and the same number of MLCPs per iteration. from_origin needs the
saved origin `a` (one vector per DS) and one product `G_j (a − q_j)` per interaction.
Both iterations converge linearly, at a rate `O(κ‖μ‖)`.

**Stopping test.** The loop continues while a constraint is violated by more than the
tolerances (§2, §3). With from_origin it also continues while
`‖q_{j+1} − q_j‖_∞ > _projectionTolq` (default 1e-12, the `fptol` of Hots).
Otherwise the loop would stop as soon as `|g| < tol`, before the normality has
converged. On the rocking block, the default `_projectionMaxIteration` is enough
(10 and 100 iterations give identical results).

**Implementation** (`ProjectionOnConstraintsTools.hpp/.cpp`, shared by both
time-stepping classes):

- `saveProjectionOrigin(nsds, a)`: before the loop, `a[ds] = q`;
- `shiftProjectionRhs(indexSet, a)`: after `computeCriteria()` (which computes
  `y(0) = g(q_j)` and `H = G(q_j)`), `y(0) += H (a − q_j)`, with `a − q_j` concatenated
  over the source and target DS of the interaction (order of the columns of `H`).
  `y(0)` is recomputed by the next `computeCriteria()`, and `MLCPProjectOnConstraints`
  is unchanged;
- `updatePositionsAfterProjection(nsds, a, mode)`: `q += p(0)` or `q = a + p(0)`.
  It returns `‖q_{j+1} − q_j‖_∞`.

Usage:

```cpp
auto s = std::make_shared<siconos::simulation::TimeSteppingCombinedProjection>(
    nsds, td, osi, osnspb_velo, osnspb_pos);
s->setProjectionUpdate(siconos::simulation::ProjectionUpdate::from_origin);
s->setProjectionTolq(1e-12);  // default
```

```python
import siconos.simulation as ss

s = ss.TimeSteppingDirectProjection(nsds, td, osi, lcp, mlcp_pos, 0)
s.setProjectionUpdate(ss.ProjectionUpdate.from_origin)
```

(`TimeSteppingCombinedProjection` has no Python binding.)

### 4.3 Metric

The projection is Euclidean (`MLCPProjectOnConstraints::_useMassNormalization = false`,
default). With `_useMassNormalization = true` the MLCP matrix is `G W⁻¹ Gᵀ` (projection
in the metric of `W`). The update, however, is still `q += Gᵀ λ` instead of
`q += W⁻¹ Gᵀ λ`, so the mass-weighted projection is currently **inconsistent**. It
would need to be fixed (update with `W⁻¹ Gᵀ λ`) before use.

---

## 5. Differences with the Scilab code Hots (and its Python port)

Hots is the Scilab code used for CMAME 2013. Its Python port, hots-python, reproduces
the Scilab results to 1e-13 and Figs. 4–7 and 9 of the article. Besides the projection
update (§4.2), the following conventions differ. Each one changes the results of the
rocking block (`RockingBlockTS-*` example data):

| | Siconos | Hots | effect on the rocking block |
|---|---|---|---|
| projection update (both algorithms) | `incremental` (default) or `from_origin` | from origin (fixed point, `fptol = 1e-12`) | 1.1e-6 at the first impact, 6e-6 over the run (Alg. 2) |
| forecast of the activation rule (Alg. 1) | linear: `g + γ h G v` | nonlinear: `g(q_k + γ h v_k)`, rule (61); rule (62) "free" forecast available | with the threshold below, the activation of a contact flips at t = 1.35 (nonlinear forecast gap 4.0e-7): max difference 1.6e-1 afterwards |
| `γ` (Alg. 1) | 1/2 | 1 (default of `projected_moreau_step`) | different impact sequence from the first impact: max difference 6.2 (the former `RockingBlockTSProj.ref` was computed with γ = 1) |
| activation threshold (Alg. 1) | 1e-7 | 1e-10 | same flip at t = 1.35 (with a threshold ≥ 1e-6, the nonlinear forecast gives the Siconos result) |
| Jacobian of the velocity Newton law (Alg. 1) | `G(q_k)` | `G` at the forecast point `q_k + γ h v_k` | first difference at the impact at t = 0.85, max difference 1.5 |
| projection (Alg. 1) | iterated until `g ≥ −tol` | one linearization (`fpitermax = 0`) by default | none on this example (one iteration was always enough) |
| activation (Alg. 2) | `g < −_constraintTolUnilateral` (1e-8 by default) | `g < −tolactive` (1e-10) | not measured (the comparisons use 1e-8 in both codes) |

With the Siconos conventions reproduced in Python, the two codes agree to round-off
(next section).

---

## 6. Validation

### 6.1 Procedure

hots-python is run on the time grid of the Siconos output, with a model built from the
data of the example. For the rocking block of `RockingBlockTS-*.cpp`:

```python
import math, numpy as np
from hots.model import Model
from hots.models import _zero_fint, _zero_kc
from hots.projected import combined_projected_moreau_fixed

LB, HB = 0.2, 0.1  # contact geometry of the example (LBlock, HBlock)
I = (1.5**2 + 1.0**2) / 12  # inertia (HeightBlock = 1.5, LengthBlock = 1)
M = np.diag([1.0, 1.0, I])
fext = np.array([0.0, -9.81, 0.0])


def g(q):
    c, s = math.cos(q[2]), math.sin(q[2])
    return np.array(
        [q[1] - 0.5 * LB * s - 0.5 * HB * c, q[1] + 0.5 * LB * s - 0.5 * HB * c]
    )


def G(q):
    c, s = math.cos(q[2]), math.sin(q[2])
    return np.array(
        [[0, 1, -0.5 * LB * c + 0.5 * HB * s], [0, 1, 0.5 * LB * c + 0.5 * HB * s]]
    )


model = Model(
    name="RB",
    q0=[0, 1.0, 0.2],
    v0=np.zeros(3),
    mass=lambda q, t: M,
    fint=_zero_fint(3),
    fext=lambda t: fext,
    nabla_fint=_zero_kc(3),
    h=g,
    nabla_h=G,
    e=0.5,
    is_linear=False,
)

dat = np.loadtxt("RockingBlockTS-Combined-FromOrigin.dat")
tr = combined_projected_moreau_fixed(
    dat[:, 0], model, theta=0.50001, tolactive=1e-8, fpitermax=100, fptol=1e-12
)
err = np.abs(dat[:, 1:7] - np.hstack([tr.q, tr.v])).max()
```

For Algorithm 1 with the Siconos conventions, the step of `hots/projected.py` is
modified accordingly: linear forecast, `γ = 1/2`, `G(q_k)` in the velocity problem,
threshold 1e-7, and a projection only when a constraint is violated by more than 1e-8.
For the `incremental` variant, the projection loop of Python is replaced by
`q ← q + G(q)ᵀλ, G Gᵀ λ = −g(q)` until `|g| < 1e-8`.

### 6.2 Results

| example (siconos-tutorials) | variant | compared with | max error on (q, v) |
|---|---|---|---|
| `BouncingBallTS-MoreauJeanCombinedProjectionOSI` (e = 0.9, h = 5e-3, T = 10) | incremental | Hots Alg. 2 | 1e-10 (q), 2e-14 (v); identical to the existing ref |
| same | from_origin | ref | 9e-16 (linear constraint) |
| `ColumnOfBeadsTS-CombinedProjection` | from_origin | ref | 9e-16 (linear constraints) |
| `RockingBlockTS-CombinedProjection` | incremental | Hots Alg. 2 with incremental update | 2.2e-15 |
| same | incremental | Hots Alg. 2 (exact projection) | 6.3e-6 |
| same | **from_origin** | **Hots Alg. 2** | **4.0e-15** |
| `RockingBlockTS-ProjectOnConstraints` (Alg. 1, level 0) | incremental | Alg. 1 with γ = 1/2, `G(q_k)`, single linearization | 2.2e-15 |
| same | **from_origin** | Alg. 1, Siconos conventions, exact projection | **2.2e-15** |

Before the fixes of §3, the bouncing ball (combined) stayed on the ground after the first
impact (error 7.5 on the ref), and the rocking block (combined) differed from Hots by
0.157.

Reference files:

| file | content |
|---|---|
| `RockingBlock/RockingBlockTS-Combined.ref` | Alg. 2, incremental |
| `RockingBlock/RockingBlockTS-Combined-FromOrigin.ref` | Alg. 2, from_origin |
| `RockingBlock/RockingBlockTSProj.ref` | Alg. 1, incremental |
| `RockingBlock/RockingBlockTSProj-FromOrigin.ref` | Alg. 1, from_origin |
| `BouncingBall/BouncingBallTS-MoreauJeanCombinedProjectionOSI.ref` | Alg. 2 (unchanged) |

The two RockingBlock examples run both variants (`simulate(ProjectionUpdate, name)`).

---

## 7. Known limitations

- `ProjectionUpdate::from_origin` is implemented for Lagrangian systems only
  (`LagrangianDS`, `LagrangianSparseDS`, Lagrangian relations). With `NewtonEulerDS`, the
  positions contain a quaternion that is normalized after each update (`normalizeq`),
  which breaks the fixed point; an exception is thrown.
- The mass-weighted projection (`_useMassNormalization`) is inconsistent (§4.3).
- The new members (`_projectionUpdate`, `_projectionTolq`) are not in the generated
  serialization lists (`io/src/generation/SiconosFull*Generated.hpp`).
- Other examples using the combined scheme still fail against their refs
  (`SliderCrankMoreauJeanCombinedProjectionOSI`, `JointsTests_*_MoreauJeanCombinedProjection`,
  NewtonEuler), as on the main branch, as does `ColumnOfBeadsTS-ProjectOnConstraints`;
  they have not been compared with a reference implementation.
