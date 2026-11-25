# Elastoplastic FEM solids: Moreau–Jean scheme with stress unknowns

This note describes how `mechanics/fem/native` integrates a linear elastic solid with
stress-based plasticity constraints (`SolidLinearTIDS` + `StressLinearTIR` +
`MohrCoulombPlasticityNSL`), the sign convention of the iteration matrix and of the
plastic variables, and how the scheme was validated. It documents the fix of
November 2026 (example `plasticBernoulliBeam` failing: the plasticity problem had no
physical solution).

- Dynamical system: `SolidLinearTIDS` (unknowns `x = [v; σ]`, `real_size() = ndof + nσ`).
- Relation: `StressLinearTIR`, `y = C σ + e`, input `Cᵀλ`.
- Law: `MohrCoulombPlasticityNSL(c, φ, 3)`: `b = c` on the first component of `y`, μ = 0
  (φ not used yet).
- Integrators: `fem::integrators::MoreauJeanGOSI` (global, used by the examples) and
  `fem::integrators::MoreauJeanOSI` (non-global; untested with plasticity).
- One-step problem: `fem::nonsmooth_formulations::GlobalFrictionContact`.

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

## 7. Open points

- The non-global `MoreauJeanOSI` path with stress interactions is untested.
- φ (Mohr–Coulomb cone) is not used.
- Energy balance of `ballAgainstPlasticBernoulliBeam`.
