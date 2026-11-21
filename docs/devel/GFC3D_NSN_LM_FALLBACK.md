# GFC3D_NSN_AC: Levenberg–Marquardt fallback and default formulation

Solver: `gfc3d_nonsmooth_Newton_AlartCurnier` (`numerics/src/FrictionContact/gfc3d_nonsmooth_Newton_AlartCurnier.c`),
solver id `SICONOS_GLOBAL_FRICTION_3D_NSN_AC`.

## 1. The method

The global friction-contact problem (velocity v, local velocity u, reaction r)

```
M v - H r = f,      u = H^T v + w,      (u, r) satisfy the Coulomb law (3x3 per contact)
```

is written as a system of nonsmooth equations Ψ(x) = 0, x = (v, u, r), where the Coulomb law of
each contact is replaced by a 3×3 function F(r, u) (`compute_ac_fun_3x3`). The semi-smooth Newton
iteration solves

```
J(x_k) d_k = -Ψ(x_k),     x_{k+1} = x_k + α_k d_k,
```

with α_k given by a Goldstein–Price line search on the merit function q(x) = ½‖Ψ(x)‖²
(`global_linesearch_goldstein_price`), or α_k = 1 without line search.

Two formulations of F are available (`SICONOS_FRICTION_3D_NSN_FORMULATION`):

- **Alart–Curnier** (`..._ALARTCURNIER_STD`): radius of the friction disk μ·max(0, r_N − ρ_N u_N)
  (projected normal reaction);
- **Jean–Moreau** (`..._JEANMOREAU_STD`): radius μ·r_N (current normal reaction), hence an extra
  coupling ∂F_T/∂r_N = −μ (r_T − ρ_T u_T)/‖r_T − ρ_T u_T‖ outside the disk.

The local Jacobians of both were checked branch by branch (normal part open/closed, tangential
part inside/outside the disk): they are the derivatives of F.

Optional shift regularization (`SICONOS_FRICTION_3D_NSN_REGULARIZATION`, `..._REGUL` = ε):
the step solves (J + εI) d = −Ψ.

## 2. The failure (FCLIB collection, `gfc3d_TEST_NSN_COLLECTION_FCLIB`)

On `Spheres-i099-356-679` (356 contacts, 12000 dofs) and `Box_Stacks-i0122-82-5` (82 contacts),
with tolerance 1e-5, ε = 1e-6 and 1000 iterations, the Jean–Moreau presets with line search did
not converge (error ~2e-4 after 1000 iterations), while every Alart–Curnier preset with line
search converged in 4–20 iterations. The failure was already present before the 2026-11 changes.

Iteration log (Jean–Moreau, line search, ε = 1e-6, Spheres): normal progress for 24 iterations,
then at iteration 24 a full step brings ‖Ψ‖ from 1.4 to 2.6e-4 and the contact statuses freeze
(90 open, 222 stick, 44 slip); afterwards ‖Ψ‖ stays at 2.134e-4 for hundreds of iterations, the
line search accepting only steps of 1e-5…2e-4 (1e-16 without regularization).

Diagnosis:

- the Jacobian is nearly singular at that point: without regularization the Newton system is
  solved only to ‖J d + Ψ‖ ≈ 1e-4…2e-6 for ‖Ψ‖ ≈ 2e-4 (LU on a nearly singular matrix);
- with the shift (J + εI) d = −Ψ, the components of d along the directions of the small
  singular values σ of J are multiplied by ~1/ε: ‖d‖ ≈ 0.1, a step that crosses many kinks of
  the nonsmooth function, so the linear model is useless and the line search only accepts tiny
  steps: the iterate is stuck at a point where q no longer decreases although Ψ ≠ 0;
- the outcome depends erratically on ε (0: stall, 1e-10: stall, 1e-8: converged in 26
  iterations, 1e-6: stall): the iterates fall or not in the basin of the solution.

A plausible origin of the near-singularity is redundant sticking contacts (each sticking contact
acts as 3 bilateral constraints u = 0; several on the same body make the reactions non unique);
not verified on these problems.

## 3. The Levenberg–Marquardt fallback

When the Newton iteration stalls, i.e. ‖Ψ_k‖ > 0.99 ‖Ψ_{k-1}‖
(`GFC3D_NSN_AC_LM_STALL_RATIO`), the direction of the iteration is computed from

```
(J^T J + λ I) d = -J^T Ψ,     λ = c ‖Ψ_k‖     (c = SICONOS_FRICTION_3D_NSN_LM_COEF)
```

(`gfc3d_nsn_ac_lm_direction`), then the usual line search is applied. Properties:

- d is a descent direction of q for any λ > 0: ∇q·d = (J^T Ψ)·d = −d^T (J^T J + λI) d < 0;
- along a singular direction of J (singular value σ) the component of d is
  σ/(σ² + λ) × (component of Ψ): damped instead of amplified (1/ε for the shift);
- λ ∝ ‖Ψ‖ vanishes at a solution, so the fast local convergence of Newton is recovered as soon
  as the iteration progresses again (the next iterations use the Newton direction).

The fallback only acts on stalled iterations: when Newton progresses, the iterates are unchanged.
Its cost is the product J^T J and its LU factorization (instead of the LU of J) on the stalled
iterations only. If the factorization fails, the Newton direction is used.

Options (`FrictionContact_options.h`):

| option | index | default | meaning |
|---|---|---|---|
| `SICONOS_FRICTION_3D_NSN_LM_FALLBACK` | iparam 23 | 1 | 0: off, 1: on |
| `SICONOS_FRICTION_3D_NSN_LM_COEF` | dparam 6 | 1.0 | c in λ = c‖Ψ‖ |

Limitation: without line search (`SICONOS_FRICTION_3D_NSN_LINESEARCH_NO`) the iteration is not
globalized and the fallback does not help (Box_Stacks, Jean–Moreau, no line search: full steps
cycle around ‖Ψ‖ ≈ 1.1e-3 with or without it).

## 4. Default formulation

The default formulation of `GFC3D_NSN_AC` is now **Alart–Curnier** (it was Jean–Moreau). This
changes the default behaviour of the solver (e.g. the python `GlobalFrictionContactProblem.solve`
helper without options); the Jean–Moreau formulation remains available with
`SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD`.

## 5. Results (macOS arm64, tolerance 1e-5, ε = 1e-6, max 1000 iterations)

`gfc3d_TEST_NSN_COLLECTION_FCLIB` (9 presets × 2 problems):

| case | problem | preset | before | after |
|---|---|---|---|---|
| 4 | Box_Stacks | JM, no line search | error 1.1e-3 (1000 it.) | unchanged, expected failure |
| 10 | Spheres | AC, no line search | error 1.9 | unchanged, expected failure |
| 12 | Spheres | JM, Goldstein–Price | error 2e-4 (1000 it.) | converged (26 it., error 1.2e-10) |
| 13 | Spheres | JM, no line search | error 4.9 | unchanged, expected failure |
| 16 | Spheres | JM, line search, warm start auto | error 2e-4 (1000 it.) | converged |
| others | | | converged | converged (unchanged) |

Run time of the collection: 39 s → 7.8 s. All the other tests (194 with FCLIB and HDF5, 181 in
the release configuration) pass.

Note: the FCLIB collections need `-DWITH_FCLIB=ON -DWITH_HDF5=ON`; `WITH_HDF5` is not declared as
an option in the CMake files (only tested), so it must be given explicitly.
