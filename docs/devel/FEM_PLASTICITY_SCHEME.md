# Elastoplastic FEM solids: Moreau–Jean scheme with stress unknowns

This note describes how `mechanics/fem/native` integrates a linear elastic solid with
stress-based plasticity constraints (`SolidLinearTIDS` + `StressLinearTIR` +
`MohrCoulombPlasticityNSL`), the sign convention of the iteration matrix and of the
plastic variables, and how the scheme was validated. It documents the fix of
November 2026 (example `plasticBernoulliBeam` failing: the plasticity problem had no
physical solution).

- Dynamical system: `SolidLinearTIDS` (unknowns `x = [v; σ]`, `real_size() = ndof + nσ`).
- Relation: `StressLinearTIR`, `y = C σ + e`, input `Cᵀλ`.
- Laws: `MohrCoulombPlasticityNSL(c, φ, 3)`: `b = c` on the first component of `y`, μ = 0
  (φ not used yet); `VonMisesPlasticityNSL(σ_y, hypothesis)` (§7);
  `DruckerPragerPlasticityNSL(η, θ, c')`, `mohrCoulombPlaneStrain(c, φ, ψ)` (§8).
- Integrators: `fem::integrators::MoreauJeanGOSI` (global, used by the examples) and
  `fem::integrators::MoreauJeanOSI` (non-global; untested with plasticity).
- One-step problems: `fem::nonsmooth_formulations::GlobalFrictionContact` (box bounds of
  `MohrCoulombPlasticityNSL`) and `fem::nonsmooth_formulations::Plasticity` (Von Mises and
  Drucker–Prager, condensation on the stresses and numerics `PlasticityProblem`, §7, §8).

---

## 1. Continuous model

```
M v̇ + Bᵀσ = F_ext + p                    (equilibrium, B: strain operator)
S σ̇ = B v − ε̇_p                           (constitutive law, S = D⁻¹: compliance)
ε̇_p = λ̇ ∂f/∂σ,  0 ≤ λ̇ ⊥ −f(σ) ≥ 0          (associated flow rule)
```

| symbol | meaning | in Siconos |
|---|---|---|
| `M` | mass matrix | `SolidLinearTIDS::mass()` |
| `B` | strain–displacement operator (element integrated) | `BMatrix()` |
| `S` | compliance, inverse of the elasticity tensor (element integrated) | `elasticityMatrix()` |
| `σ` | generalised stresses (6 per beam element) | `stress()` |
| `ε_p`, `ε̇_p` | plastic deformation and its rate | `plasticDeformation()`, `plasticRate()` |

A plasticity interaction is `y = C σ + c ≥ 0` (yield function `f = −y`). The example
uses two interactions per element with `C = +e₁ᵀ` and `C = −e₁ᵀ` to get `|σ₁| ≤ c`.

## 2. Time discretisation (θ-method)

With `x_{k+1} = [v_{k+1}; σ_{k+1}]`:

```
W x_{k+1} = q + H r,      W = [  M      hθBᵀ ]
                              [ −hθB     S   ]

q = [ M v_k − h(1−θ) Bᵀσ_k + h F_ext,θ ]      (MoreauJeanGOSI::computeResidu)
    [ S σ_k + h(1−θ) B v_k             ]

H r = [ p ; Cᵀλ ],   Cᵀλ = −h ε̇_p            (stress rows: plastic increment)
y   = Hᵀ x + b,      b = c                     (GlobalFrictionContact)
```

The second row is the discretised constitutive law
`S σ_{k+1} − hθ B v_{k+1} = S σ_k + h(1−θ) B v_k − h ε̇_p`.

### Why `[−hθB, S]` and not `[hθB, −S]`

Both write the same equations; the second one (used until November 2026) makes `W`
symmetric but **indefinite**. Its inverse has a negative definite stress block
`−(S + h²θ² B M⁻¹ Bᵀ)⁻¹`, so the Delassus operator of a stress interaction
`Hᵀ W⁻¹ H` is negative: increasing λ moves `y` the wrong way. On
`plasticBernoulliBeam` at the first plastic step:

```
W_loc (normal block) = [[−2.72e10,  2.72e10], [ 2.72e10, −2.72e10]],  eig ≤ 0
q_loc = [4.06e9, −6.3e7]
```

The only LCP solution activated the "+" interaction (σ ≥ −c) and sent σ to −c; the
solvers (NSGS_WR, NSN_AC, VI_EG, IPM) failed or diverged, and the stress went through the
yield bound. With `[−hθB, S]` the symmetric part of `W` is `diag(M, S)` (positive
definite), the Delassus operator is positive and the multiplier pushes σ back to the
yield surface. `W` is not symmetric, which the global solvers (LU of `M` in NSGS_WR)
handle.

The non-global integrator uses the same `W` in increment form
(`W Δx = rhs`, `Δx = x_{k+1} − x_k`):

```
rhs = [ −h Bᵀσ_k + F_ext,θ + p ;  h B v_k − h ε̇_p ]
```

(`MoreauJeanOSI::computeResidu`, `integrate`, `computeIteration`).

## 3. Plastic variables

- `StressLinearTIR::computeInput` **accumulates** `Cᵀλ` into `plasticRate()` (the "+"
  and "−" interactions share the same block), like `p` for contacts.
- `SolidLinearTIDS::resetAllNonSmoothParts` / `resetNonSmoothPart(1)` zero it before
  each input update (`Simulation::updateAllInput`). Without this reset the "rate" was
  the running sum of all multipliers (it never went back to 0 on unloading).
- `fem::integrators::MoreauJeanOSI::updateInput` converts it, after the relations, into
  the physical rate `ε̇_p = −Cᵀλ / h` (positive when σ is on the +c surface).
- `MoreauJeanOSI::updateState` integrates `ε_p,k+1 = ε_p,k + h ε̇_p` from
  `plasticDeformationMemory()` (idempotent over Newton iterations); `swapInMemory`
  stores `ε_p`.

## 4. Friction coefficient

`GlobalFrictionContact::updateMu` (and `compute_nslaw_contribution`) push μ = 0 for a
`MohrCoulombPlasticityNSL`: the law only bounds the first component of `y` by `c`.
Before, nothing was pushed and the solver read an empty `_mu` vector (uninitialised
memory). A Mohr–Coulomb cone would use μ = tan φ together with the cone `b`, not done.

## 5. Discrete operators of continuum elements (T3)

For a 2D continuum element e of area A_e and thickness t (`Material` thickness):

- `BMatrix()`: the strain–displacement matrix scaled by `sqrt(A_e t)`, so that the stress
  unknowns are `σ̃_e = sqrt(A_e t) σ_e` (σ_e the physical stress, Voigt order
  (σxx, σyy, σxy), engineering shear strain).
- `elasticityMatrix()`: block diagonal, `S_e = D⁻¹` with D the plane strain or plane stress
  tensor of the material (`elasticityTensor2D`). With this scaling, `K = Bᵀ S⁻¹ B` (checked in
  `SolidLinearTIDSTest::testPlaneStressCompliance`). A yield criterion `f(σ) ≤ 0` written on
  the unknowns becomes `f(σ̃_e / sqrt(A_e t)) ≤ 0`; for a positively homogeneous criterion (von
  Mises, Drucker–Prager without cohesion) the yield stress of element e is
  `sqrt(A_e t) σ_y`.
- Mass: integral of `ρ t Nᵀ N` (consistent, default) or its row-sum lumping
  (`MassType::lumped`, constructor argument of `SolidLinearTIDS`). The lumped mass is
  diagonal, which keeps the stress Schur complement `S + h²θ² B M⁻¹ Bᵀ` sparse.

The analysis type is a constructor argument of `Material`:
`Material(ρ, E, ν, AnalysisType2D::plane_stress, t)` (the 3-argument constructor keeps plane
strain, t = 1). Before November 2026 the compliance was always the plane strain one (the
stiffness already used the analysis type) and the thickness did not enter the mass.

These operators are identical to those of the pytripopfem code (native T3 backend, one Gauss
point, `scaling_sigma_vector = sqrt(weight det J)`, lumped mass) to round-off on the strip with
a circular hole (`strip_with_circular_hole_72.msh`, plane stress): relative differences
6.7e-16 (B), 1.4e-16 (S), 4.5e-16 (M).

## 6. Validation (macOS arm64, release)

`plasticBernoulliBeam` (h = 1e−5, 1000 steps, c = 2e9):

| | before | after |
|---|---|---|
| solver warnings | 213 | 0 |
| max σ₁₇ | 2.48e9 (> c) | 2.000000000000e9 |
| plastic phase | none | steps 265–405, ε̇_p > 0 only on the yield surface |
| ε̇_p during unloading | – | exactly 0 |
| final ε_p | 0 | 0.1745 |

The stresses agree to ~3 digits with the June 2026 reference (computed with an older
NSGS_WR that happened to return the clamped solution); its "rate" column was −Σλ. The
reference has been regenerated and the example compares stresses and plastic variables
with tolerances relative to their magnitude.

`ballAgainstPlasticBernoulliBeam`: 602 → 0 solver warnings, |σ| ≤ c, plastic dissipation
≥ 0 (its energy post-processing now uses `plasticDeformation`, not `plasticRate`). The
total energy still drifts (also before the change): open.

Non-regression: `bernoulliBeam`, `ballAgainstBernoulliBeam`, `TH4`, `T3`, `T3-GOSI`,
`T3_collision_detection(_with_friction)`, `Hertz_T3_collision_detection` give
bit-identical results (without plasticity, negating a row and its right-hand side does
not change the solution).

## 7. Von Mises plasticity: condensed one-step problem

### Law and relation

`VonMisesPlasticityNSL(σ_y, hypothesis)` (kernel): the output y of a size-3 interaction is a
stress in Voigt notation (σxx, σyy, σxy) that must stay in the elastic domain K, and the
multiplier is in the opposite of the normal cone (associated flow rule):

    0 ∈ λ + N_K(y),   Cᵀλ = −h ε̇_p

- plane stress: K = {σxx² + σyy² − σxx σyy + 3 σxy² ≤ σ_y²};
- plane strain (incompressible plastic flow, σzz = (σxx + σyy)/2):
  K = {(σxx − σyy)²/4 + σxy² ≤ σ_y²/3}.

One interaction per element, relation `StressLinearTIR(C_e)` with
`C_e = SolidLinearTIDS::elementStressSelection(e)` (the 3 stresses of the element divided by
`sqrt(A_e t)`, §5): y is then the physical stress and one law with the physical σ_y serves
all the elements. `b = 0` (no contribution of the law to the output).

### Condensation (`fem::nonsmooth_formulations::Plasticity`)

The friction solvers cannot be used: in `GlobalFrictionContact` the cone constrains the
multiplier (and the De Saxcé term modifies the condition on y), whereas plasticity constrains
the output y. The problem assembled by `GlobalFrictionContact::preCompute` (W of
`MoreauJeanGOSI`, boundary conditions included),

    W x = q + H z,   y = Hᵀ x + b,   0 ∈ z + N_K(y),   x = [v; σ],

is condensed on the constrained unknowns c (rows of x touched by H; r the others: velocities,
unconstrained stresses). H_c (H restricted to c) must be square and invertible (each
constrained stress belongs to exactly one interaction):

    W_red = W_cc − W_cr W_rr⁻¹ W_rc,   q̃ = q_c − W_cr W_rr⁻¹ q_r,   W_red x_c = q̃ + H_c z,

and with `x_c = H_c⁻ᵀ (y − b)`:

    z = M_p y + q_p,   M_p = H_c⁻¹ W_red H_c⁻ᵀ,   q_p = −M_p b − H_c⁻¹ q̃.

This is the numerics `PlasticityProblem` `0 ∈ M_p y + q_p + N_K(y)` (Von Mises model,
`PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT`, hypothesis of the laws), solved by
`plasticity_2d_driver` (default `PLASTICITY_2D_NSGS`). Then `x_r = W_rr⁻¹(q_r − W_rc x_c)`,
`_globalVelocities = x`, `_z = z` (λ of the interactions), `_w = y` (warm start of the next
step). When W_rr is diagonal (lumped mass and every stress constrained: W_rr = M with the
boundary rows) it is inverted entry-wise and W_red = S + h²θ² B M⁻¹ Bᵀ (restricted) is
sparse; otherwise a sparse LU of W_rr is used (the result W_rr⁻¹ W_rc is dense in general).
For the lumped T3 case, M_p is pytripopfem's W (up to the scaling of the stresses).

W and H are time-invariant (linear elasticity, `StressLinearTIR`): the full assembly
(`GlobalFrictionContact::preCompute`) and the condensation (selections, W_rr⁻¹ W_rc, M_p,
H_c⁻¹, the numerics matrix and the model parameters) are done at the first step by
`Plasticity::assemble` and kept in a `PlasticityCondensation`; `Plasticity::preCompute`
afterwards only recomputes q and b, and `solve_condensed` q̃, q_p and the solution. The
condensation is rebuilt when the number of interactions changes.

`Plasticity::setProblemDump(prefix, steps)` writes the numerics problem passed to the driver
(numerics `.dat` format, with the coordinates and the hypothesis) into `<prefix>_<n>.dat` for
the chosen calls n of `compute` (time steps). The problems of the strip below (steps 1, 4, 10,
100) are the numerics test data `plasticity_2d_vonmises_strip_<n>.dat`
(`plasticity_2d_strip_test`, see `docs/devel/PLASTICITY.md`).

### Prescribed velocities in `fem::integrators::MoreauJeanGOSI`

The global integrator solves for x_{k+1} itself. The right-hand side is now computed from the
unmodified blocks,

    q_v = M v_k − h(1−θ) Bᵀσ_k + h F_ext,θ,   q_σ = S σ_k + h(1−θ) B v_k,

and a prescribed velocity v_j = g_j(t_{k+1}) is eliminated with the original column of W
(stored in `iterationMatrixBoundaryConditions` before its row and column are zeroed):
`q −= g_j W(:, j)` in the velocity **and stress** rows (−hθB(:, j) g_j), then
`q_j = W_jj g_j`. Before November 2026 the increment formula of `MoreauJeanOSI`
(Δv = g − v_k) was used and the stress rows were not corrected: a non-zero prescribed
velocity was added at every step (velocity growing linearly) with the wrong sign, and the
stress equation lagged one step. Zero prescribed velocities (all the existing examples) are
not affected (results equal up to round-off: `ballAgainstBernoulliBeam` 2e-10 on a column of
order 1e5, its comparison is now relative).

### H of stress interactions with boundary conditions (`OSNSMatrix::fillH`)

`fillH` zeroed the columns of an interaction block at the constrained **velocity** indices.
For a `StressLinearTIR` the block acts on the stresses: unrelated stress components were
removed from H (15 on the strip below, element 0 of a clamped beam). Stress interactions
are now skipped.

### Validation: strip with a circular hole (pytripopfem)

Quarter of a strip (x ∈ [0, 0.5], y ∈ [0, 1], hole of radius 0.25), mesh
`strip_with_circular_hole_72.msh` (41 nodes, 60 T3), plane stress, E = 210e3, ν = 0.3,
ρ = 7.8e-9, σ_y = 210, lumped mass, symmetry conditions on the bottom and left edges, top edge
pulled at v = 100, h = 1e-6, 100 steps, θ = 1. Reference: pytripopfem native T3 backend
(`examples/Strip_with_hole`, ipopt tolerance 1e-10, fenicsx-siconos-0.10 container).

| NSGS tolerance | sweeps (total / max per step) | max rel. diff u | σ | top reaction R |
|---|---|---|---|---|
| 1e-12 | ≈ 1e7 / 1e5 (cap reached at steps 1–10) | 4.0e-6 | 2.2e-4 | 4.4e-6 |
| 1e-8 | 5.4e6 / 7.3e4 | 2.8e-5 | 3.5e-4 | 7.1e-6 |

The reaction of the top edge rises linearly (7.997 at step 1) and saturates at 56.511
(pytripopfem 56.511; net-section limit load σ_y · 0.25 = 52.5, overestimated by the coarse
mesh); 23 elements are plastic at the end in both codes; ‖σ‖_vm ≤ σ_y to 4e-16. The ipopt
solutions of the reference have a natural-map error of 1e-7–2e-6 (Siconos NSGS: the
tolerance), which bounds the stress agreement. The step problems of pytripopfem, re-solved
by `plasticity_2d_driver` (M = W, q = −q_σ, per-element σ_y scaled by sqrt(A_e)), give the
same stresses (relative difference ≤ 1e-4, ipopt accuracy).

Tutorial: `siconos-tutorials/examples/mechanics/FEM/Plasticity/StripWithHole`
(`strip_with_hole.cpp`, NSGS tolerance 1e-8; 20 steps in the tests, about 25 s on an M-series
Mac; 100 steps, about 2 min, when run directly).

Unit test `SolidLinearTIDSTest::testVonMisesTension`: plane stress uniaxial tension of
`square_6.msh`, σxx = E ε then σ_y in every element (exact for T3), σyy = σxy = 0, to 2e-7
outside the dynamic transients (start of the loading, onset of yield: 6e-4, decaying in
3 steps); von Mises stress ≤ σ_y.

### Cost

The NSGS of numerics needs 4e4–1e5 sweeps per step on the strip (about 20 µs per sweep for
60 elements). The cost is in the solver, not in the condensation: NSGS needs as many sweeps on
pytripopfem's W directly. In this quasi-static regime (h = 1e-6 ≫ wave transit time 2e-8 s
through an element) W_red is dominated by the rank-deficient term h² B M⁻¹ Bᵀ (≈ 1e-2 against
S ≈ 5e-6), condition number ≳ 1e4, and Gauss–Seidel needs O(κ) sweeps. The NSGS relaxation
option does not help: it relaxes after the exact local solve, σ ← ω σ_loc + (1 − ω) σ_old,
which leaves the (strictly convex) elastic domain for ω > 1 (no convergence for ω ≥ 1.8);
ω = 1.5 gives no gain. A tolerance of 1e-8 instead of 1e-12 saves 40 % (reaction within
1.3e-6).

The global semi-smooth Newton of numerics (`PLASTICITY_2D_NSN`, see `docs/devel/PLASTICITY.md`)
is the default solver of `fem::nonsmooth_formulations::Plasticity`. With tolerance 1e-10 it runs
the 100 steps of the strip in 0.33 s (NSGS: ~100 s), with 4–9 iterations per plastic step and a
better agreement with pytripopfem (displacements 4.0e-6, reaction 4.5e-6, stresses 2.2e-4 —
ipopt's accuracy). Finer meshes of pytripopfem (same data):

| mesh (T3 elements, stresses) | NSN iterations per plastic step (mean / max) | time | memory | final reaction |
|---|---|---|---|---|
| 60, 180 | 8 / 9 | 0.33 s | | 56.51 |
| 636, 1 908 | 13 / 18 | 17 s | 0.25 GB | 53.33 |
| 1 446, 4 338 | 19 / 30 | 104 s | 1.3 GB | 52.94 |

The reaction converges to the limit load σ_y · 0.25 = 52.5 from above.

Before 2026-11-26 the 1446-element mesh took 167 s (636: 24 s), ~70 s of them in the kernel:
- the dense LU (`FullPivLU`) of the (n_v + n_σ)² iteration matrix in
  `fem::MoreauJeanOSI::initializeIterationMatrix` (~30 s), useless for
  `fem::MoreauJeanGOSI` (the one-step problem holds W): skipped there
  (`needsIterationMatrixLU()`);
- the renumbering of the graphs (`Topology::setProperties` at each step through
  `Simulation::initialize`, `updateIndexSets` at each Newton iteration): O(V + E), and all
  the element interactions share the solid, so the interaction graphs have n_el²/2 edges.
  `SiconosGraph::update_{vertices,edges}_indices` now return at once when the graph has not
  changed since the last renumbering (flags set by its mutators);
- `OSNSMatrix::fillW`/`fillH` (scan of the dense W and of the dense C of each
  `StressLinearTIR`) and the condensation at each step: done once (see above).

The initialization of the 1446-element mesh takes 1 s, and the kernel ~3 % of a step: the
rest is the numerics NSN, 90 % of it in the sparse LU of the Newton matrix (CSparse `cs_lu`,
no supernodes). The memory is still quadratic in the number of elements (1.1 GB at
initialization for 1446 elements, 3.5 GB for 2612): dense iteration matrix (n_v + n_σ)²,
dense 3 × n_σ C of each `StressLinearTIR`, and the n_el²/2 edges of the interaction graphs.

## 8. Drucker–Prager (Mohr–Coulomb) plasticity: footing

### Law

`DruckerPragerPlasticityNSL(η, θ, c')` (kernel): y = (σxx, σyy, σxy) (tension positive) in the
cone of the numerics `PlasticityProblem` in Voigt coordinates,

    ‖σ_T‖ ≤ η σ_N + c',   σ_N = −(σxx + σyy)/√2,   σ_T = ((σxx − σyy)/√2, √2 σxy),

i.e. R ≤ −η p + c'/√2 with p = (σxx + σyy)/2 and R = √(((σxx − σyy)/2)² + σxy²) (radius of the
Mohr circle), and the flow rule −(w + θ ‖w_T‖ e_N) ∈ N_K(y) (θ = 0: associated).
`mohrCoulombPlaneStrain(c, φ, ψ)` gives the Mohr–Coulomb criterion of the plane (exact in
plane strain when σzz is the intermediate principal stress), R ≤ −p sin φ + c cos φ:

    η = sin φ,   c' = √2 c cos φ,   θ = sin φ − sin ψ   (ψ: dilatancy angle).

(The dilatancy: on the yield surface the plastic strain rate ε̇_p = −w has
ε̇_N = (θ − η) ‖ε̇_T‖ = −sin ψ ‖ε̇_T‖.) The condensed problem `fem::nonsmooth_formulations::
Plasticity` passes the Drucker–Prager model (η, θ, c' per interaction) to the numerics; the
driver removes the cohesion by a translation and maps the Voigt coordinates to the cone
coordinates. pytripopfem writes the same criterion c₀ cos φ − (sin φ/2) I − J ≥ 0 (I = σxx + σyy,
J = R): its η is sin φ / 2 and its dilatancy sin ψ / 2.

### Validation: half footing (pytripopfem `examples/Footing/half_footing.py`)

Soil [0, 2]², 639 T3 (`half_footing_2x2_b02_size01.msh`), plane strain, E = 1e7, ν = 0.3,
ρ = 1.5e3, lumped mass, c = 2e3, φ = 25°, smooth footing of half-width 0.2 pushed at
v_y = −0.01, h = 0.05, 40 steps. Reference: pytripopfem native T3 backend, its semi-smooth
Newton solver (`compute_sigma_RMA_DP_NEW_SCALED`, the NSGS warm start through the old
`sn.MohrCoulomb2DProblem` API skipped), fenicsx-siconos-0.10 container.

| ψ | max rel. diff. u | σ | footing reaction | reaction at step 40 |
|---|---|---|---|---|
| 25° (associated, θ = 0) | 3.9e-12 | 3.3e-11 | 2.1e-13 | −11 897.25 |
| 10° (θ = 0.249) | 6.9e-12 | 2.8e-11 | 4.3e-13 | −11 562.82 |

Both codes solve the step problems to ~1e-10, hence the round-off agreement; the parameter
mapping (η, c', θ) is thus validated for associated and non-associated flow. 590 elements are
on the yield surface at the end (yield function ≤ 1.5e-14 c cos φ); Prandtl's weightless
bearing capacity (−8 288 for the half footing) is overestimated by the T3 mesh and the small
domain. Cost: 13 s for 40 steps (NSN 7–13 iterations per step); NSGS reaches an error of 3e-4
after 2e4 sweeps on one of these problems. Tutorial:
`siconos-tutorials/examples/mechanics/FEM/Plasticity/Footing`; numerics test data
`plasticity_2d_druckerprager_footing_fem_40.dat`, `..._psi10_40.dat`.

## 9. Open points

- The non-global `MoreauJeanOSI` path with stress interactions is untested.
- φ (Mohr–Coulomb cone) is not used.
- Energy balance of `ballAgainstPlasticBernoulliBeam`.
- Kernel costs for large meshes (§7, Cost): dense iteration matrix of SolidLinearTIDS and
  its dense LU, assembly of W and H at each step from dense matrices, dense C of
  `StressLinearTIR`; a sparse iteration matrix and a sparse/offset relation are needed.
- `PLASTICITY_2D_NSN`: new LU of H at each iteration (reuse the symbolic factorization).
- NSGS relaxation of the Plasticity solvers: relax before the projection (projected SOR) or
  forbid ω > 1.
- `StressLinearTIR` takes a dense C of size 3 × nσ per element (memory O(n_el²)); a sparse
  or offset-based relation is needed for large meshes.
- `MoreauJeanOSI::computeIterationMatrixBoundaryConditions` prints a warning (W not
  symmetric) for every SolidLinearTIDS with boundary conditions; the elimination uses
  columns only and is valid for a non-symmetric W.
