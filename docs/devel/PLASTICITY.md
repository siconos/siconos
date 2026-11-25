# Plasticity component (`numerics/src/Plasticity`)

2D plasticity problems solved element by element (one "cone" per element): Drucker–Prager
(Mohr–Coulomb) with dilatancy and cohesion, Von Mises (J2), stresses in cone or Cartesian
coordinates.

## Files

| file | content |
|---|---|
| `PlasticityProblem.{h,c}` | problem, models (`Plasticity_DruckerPrager_model`: eta, theta, cohesion; `Plasticity_VonMises_model`: sigma_y), `stress_coordinates`, cohesion shift |
| `plasticity_2d_driver.c` | driver: Cartesian → cone change of variables, cohesion translation, trivial case, registry |
| `plasticity_2d_stress_coordinates.{h,c}` | Cartesian ↔ cone maps of stresses, strain rates, M and problems |
| `plasticity_2d_compute_error.{h,c}` | natural-map error (both models, cohesion, coordinates), `plasticity_2d_DruckerPrager_zero_is_solution` |
| `plasticity_2d_druckerprager.{h,c}` | Drucker–Prager projection local solvers |
| `plasticity_2d_onecone_nonsmooth_Newton_solvers.c` | Drucker–Prager semi-smooth Newton local solvers |
| `plasticity_2d_vonmises.{h,c}` | Von Mises domain, projection, error, return mapping, fixed point and semi-smooth Newton local solvers, their SolverOptions and registration |
| `plasticity_2d_nsgs.c`, `plasticity_2d_nsgs_generic_adapter.c` | NSGS (and generic NSGS) |
| `numerics/python/plasticity_pybind11.cpp` | Python: problems (`PlasticityProblem_new`, `PlasticityProblem_new_VonMises`), enums, `plasticity_2d_driver`, `plasticity_2d_compute_error` |
| `Plasticity_options.h` | solver ids and options |
| `test/plasticity_2d_{simple_solve,vonmises,cohesion,cartesian}_test.c` | tests (see below) |

MixedProblems uses the same models per element (see `MIXED_PROBLEMS_NSGS.md`).

## Problem and cone coordinates

For each element, the solvers work in cone coordinates sigma = (sigma_N, sigma_T1, sigma_T2).
These are the components of the 2x2 symmetric stress tensor in an orthonormal basis
(Frobenius inner product) whose first vector is along the identity:

    sigma_N  = −(sxx + syy)/sqrt2      (compression positive)
    sigma_T1 =  (sxx − syy)/sqrt2
    sigma_T2 =  sqrt2 sxy

so that:

* ||sigma_T|| = ||dev sigma|| (the in-plane deviator) = sqrt2 × (radius of the Mohr circle);
* sigma_N = −sqrt2 × (in-plane mean stress).

In mechanics terms, cone coordinates are **hydrostatic–deviatoric coordinates**:

* sigma_N is the hydrostatic (spherical) part of the stress;
* (sigma_T1, sigma_T2) are the components of the deviator in an orthonormal basis;
* the basis is the 2D Cartesian form of the Haigh–Westergaard frame. The Haigh–Westergaard
  (or Lode) coordinates (xi, rho, theta) are its cylindrical version.

The code keeps the name "cone" (`PLASTICITY_STRESS_COORDINATES_CONE`, `*_to_cone`) and the
N/T indices. They come from the second-order cone of the Drucker–Prager model and the friction
notation shared with the solvers. The coordinates themselves do not depend on the model: the
Von Mises yield set is a cylinder.


With w = M sigma + q (the strain rate, dual of the stress), the problem on each element is the
variational inequality

    sigma in K,   −(w + theta ||w_T|| e_N) in N_K(sigma),

where K is the yield set of the element and theta = 0 for Von Mises.

## Cartesian stress coordinates

### Change of variables

`PlasticityProblem::stress_coordinates`:

| value | stress s | strain rate w (dual) | P (sigma_c = P s) |
|---|---|---|---|
| `PLASTICITY_STRESS_COORDINATES_CONE` (default) | sigma_c | w_c | I |
| `..._CARTESIAN_MANDEL` | (sxx, syy, sqrt2 sxy) | (exx, eyy, sqrt2 exy) | P_M, orthogonal |
| `..._CARTESIAN_VOIGT` | (sxx, syy, sxy) | (exx, eyy, 2 exy) | P_M diag(1, 1, sqrt2) |

with

    P_M = [ −1/sqrt2  −1/sqrt2  0 ]
          [  1/sqrt2  −1/sqrt2  0 ]
          [  0         0        1 ]

The strain rate is mapped by the dual map, w_c = P^{-T} w_s, so that the power is preserved:
s·w_s = sigma_c·w_c (= sxx exx + syy eyy + 2 sxy exy).

Since N_{P^{-1}K}(s) = P^T N_K(P s), the variational inequality is invariant under this change
of variables: s solves the Cartesian problem (M_s, q_s) if and only if sigma_c = P s solves the
cone problem

    M_c = P^{-T} M_s P^{-1},   q_c = P^{-T} q_s,

with the same model parameters (eta, theta, cohesion, sigma_y) in cone coordinates.

* Mandel's notation: P^{-T} = P, so M_c = P M_s P^T has the conditioning of M_s.
* Voigt's notation: P is not orthogonal (the shear row and column of each block are scaled by
  1/sqrt2 in M_c). The change of variables is still exact.

### Implementation: once and for all in the driver

The driver performs a general change of variables s = T y on every element
(`plasticity_2d_change_of_variables`), y being the unknown of the solved problem:

    T = P^{-1} D,   M_y = T^T M T,   q_y = T^T q,   w_y = T^T w,

with P the Cartesian → cone map (identity in cone coordinates) and D the plane stress
Von Mises scaling (identity otherwise, see "Plane stress: the ball model" below). Both are
composed into one 3x3 map, so M is transformed only once.

When T != I, `plasticity_2d_driver`:

1. builds M_y = T^T M T once (`plasticity_2d_matrix_change_of_variables`):
   * dense M: 3x3 blocks transformed in place in a copy;
   * sparse-block M with block sizes multiple of 3: each block of a copy is transformed, and
     the structure is kept;
   * otherwise (NM_SPARSE): products with sparse block-diagonal matrices, giving an NM_SPARSE
     result;
2. computes q_y = T^T q, and y = T^{-1} s for the initial guess;
3. solves the problem in y, a shallow copy of the problem that shares the model parameters
   (with the plane stress hypothesis replaced by the ball), so the cohesion translation and
   the trivial-case tests apply unchanged;
4. maps the solution back: s = T y, w = T^{-T} w_y.

The other entry points:

* The solvers and local solvers are not changed.
* `plasticity_2d_nsgs` and `plasticity_2d_nsgs_generic` called directly on a Cartesian or
  plane stress problem return an error (`info = 1`).
* `plasticity_2d_compute_error` accepts both coordinates. It computes w = M s + q in the given
  coordinates, then the natural-map error element by element on (P s, P^{-T} w).
* `plasticity_2d_problem_to_cone` returns a deep copy in cone coordinates;
  `plasticity_2d_problem_change_of_variables` returns a deep copy of the problem in y (the one
  the driver solves).
* `plasticity_printInFile` does not save the coordinates.

### Mechanical meaning in Cartesian coordinates

* **Drucker–Prager / Mohr–Coulomb**: with eta = sin(phi), ||sigma_T|| ≤ eta sigma_N + c is the
  in-plane Mohr–Coulomb criterion R_Mohr ≤ −sin(phi) p + c/sqrt2 (p the in-plane mean stress,
  tension positive). It is exact in plane strain when szz is the intermediate principal
  stress.
* **Von Mises**: ||sigma_T|| ≤ sqrt(2/3) sigma_y is the J2 criterion with
  szz = (sxx + syy)/2, i.e. plane strain with incompressible (rigid-plastic) flow:
  * pure shear yields at sxy = sigma_y/sqrt3;
  * plane strain uniaxial loading yields at sxx − syy = 2 sigma_y/sqrt3.
  This is the default hypothesis, `PLASTICITY_VONMISES_PLANE_STRAIN`.
* **Von Mises, plane stress** (szz = 0, `PLASTICITY_VONMISES_PLANE_STRESS`): see "Plane
  stress: the ball model" below.
* **Not covered**: elasto-plastic plane strain with nu < 1/2, where szz is an independent
  unknown (4 stress components).

## Drucker–Prager (Mohr–Coulomb)

### Model

    K = { ||sigma_T|| <= eta sigma_N + c },

* eta: cone coefficient;
* theta: dilatancy coefficient, entering through the term theta ||w_T|| e_N;
* c >= 0: the cohesion (`Plasticity_DruckerPrager_model::cohesion`), one value per cone, owned
  by the problem. NULL means no cohesion.

### Trivial solution

sigma = 0 solves a cone without cohesion if and only if q + theta ||q_T|| e_N is in the dual
cone:

    eta ||q_T|| <= q_N + theta ||q_T||     (plasticity_2d_DruckerPrager_zero_is_solution)

q_N >= 0 is enough only when theta >= eta. Both the driver (`plasticity_2d_checkTrivialCase`)
and the Newton local solver use this test.

### Cohesion: translation of the problem

With a = c/eta (eta > 0 is required when c != 0), K = K_eta − a e_N, where K_eta is the cone
without cohesion. For sigma~ = sigma + a (a = a e_N on each cone):

    w = M sigma + q = M sigma~ + q~,   q~ = q − M a,

and N_K(sigma) = N_{K_eta}(sigma~): sigma~ solves the problem without cohesion with q~.

* `plasticity_DruckerPrager_cohesion_shift` builds a. It returns NULL if there is no non-zero
  cohesion, and an error if eta <= 0 with c != 0.
* `plasticity_2d_driver`:
  * translates the problem (q~, and sigma~ = sigma + a as initial guess);
  * hides the cohesion during the solve, so the solvers and their error computation see the
    problem without cohesion;
  * solves, then translates back (sigma = sigma~ − a; q and the cohesion are restored).

  A solver called directly (e.g. `plasticity_2d_nsgs`) ignores the cohesion: use the driver.
* `plasticity_2d_compute_error` takes the cohesion into account (error on sigma + a).
* MixedProblems (`mixedProblem_set_druckerprager_cohesion`): the translation is done for each
  element in the element-wise Gauss–Seidel (q_e~ = q_e − M_ee a e_N, z_e~ = z_e + a e_N), and
  in `mixedProblem_compute_error`.
* Python: `PlasticityProblem_new(dim, M, q, eta, theta, cohesion=None)`.
* The data files do not store the cohesion.

## Von Mises (J2)

### Yield domain

`Plasticity_VonMises_model::hypothesis` (`PlasticityVonMisesHypothesis`) selects the yield set:

| hypothesis | yield set in cone coordinates | meaning |
|---|---|---|
| `PLASTICITY_VONMISES_PLANE_STRAIN` (default) | cylinder ‖sigma_T‖ ≤ sqrt(2/3) sigma_y | plane strain, szz = (sxx + syy)/2 |
| `PLASTICITY_VONMISES_PLANE_STRESS` | ellipsoid sigma_N²/2 + (3/2)‖sigma_T‖² ≤ sigma_y² | plane stress, szz = 0 |
| `PLASTICITY_VONMISES_BALL` | ball ‖sigma‖ ≤ sigma_y | plane stress in the scaled variables |

The rest of this section first describes the plane strain cylinder. The plane stress model is
described after the iterative solvers.

With the 3D J2 convention sigma_eq = sqrt(3/2) ||dev sigma||, the plane strain Von Mises domain
is the cylinder

    C = { sigma : ||sigma_T|| <= R },   R = sqrt(2/3) sigma_y,

where `sigma_y` (parameter of `Plasticity_VonMises_model`) is the **uniaxial yield stress**.

The previous implementation had three problems:

* it used Cartesian components (sigma_x, sigma_y, tau_xy) with the plane-stress formula;
* it applied a radial scaling of the whole stress, which is not a projection;
* the data file `plasticity_2d_vonmises_example1.dat` used Drucker–Prager-like data.

### Problem

The normal cone N_C(sigma) is:

* {0} inside C;
* {lambda sigma_T/R, lambda >= 0} on the boundary, with no N component.

Hence, with sigma in C and −w in N_C(sigma), equivalently sigma = P_C(sigma − rho w) for any
rho > 0:

* w_N = 0: the plastic flow is isochoric, and sigma_N is not constrained by the yield
  condition;
* w_T = −lambda sigma_T / R, lambda >= 0, lambda (R − ||sigma_T||) = 0.

The Euclidean projection P_C is the radial return of sigma_T, with sigma_N unchanged
(`plasticity_2d_projectionOnVonMises`). The natural-map residual sigma − P_C(sigma − w) is used
by `plasticity_2d_compute_error` (absolute error, as for Drucker–Prager) and by
`mixedProblem_compute_error`.

sigma = 0 lies inside C, so it is a solution only if q = 0.

### Return mapping: exact one-element solver (`plasticity_2d_vonMises_return_mapping`)

When M is symmetric positive definite, the one-element problem is the optimality condition of
min 1/2 sigma^T M sigma + q^T sigma over C. Its solution is the closest-point projection, in
the M-norm (M is a compliance), of the trial stress −M^{-1} q on C, i.e. the classical return
mapping. The solver computes it exactly, also for a non-symmetric M.

For a 3x3 local matrix M (column-major), w_N = 0 gives

    sigma_N = −(q_N + M_NT sigma_T) / M_NN

and a 2D problem on the disk D_R:

    S sigma_T + q~ in −N_{D_R}(sigma_T),   S = M_TT − M_TN M_NT / M_NN,
                                            q~ = q_T − M_TN q_N / M_NN.

* Elastic case: sigma_T = −S^{-1} q~, if its norm is <= R.
* Plastic case: (S + mu I) sigma_T = −q~ with ||sigma_T|| = R and mu = lambda/R > 0.
  * The scalar equation g(mu) = 1/||sigma_T(mu)|| − 1/R = 0 is increasing, and concave for a
    symmetric S (as in trust-region subproblems).
  * It is solved by Newton's method from mu = 0, safeguarded by bisection on a bracket
    [lo, hi].
  * hi = ||q~||/R is enough when the symmetric part of S is positive semi-definite. Otherwise
    hi is doubled until ||sigma_T(hi)|| <= R.
  * Finally sigma_T is scaled to ||sigma_T|| = R.

The solver requires M_NN > 0 and S + mu I invertible. It returns 1 otherwise, and 2 if the
scalar equation is solved only to ~1e3 eps.

### Iterative one-element solvers on the natural map

These solvers use the same options as the friction local solvers (`Plasticity_options.h`,
`PLASTICITY_NSN_*`). Their stopping test is ||sigma − P_C(sigma − rho0 w)|| <= tol (1 + ||sigma||),
with rho0 = 1/||M||_2 (power iteration on M^T M).

**Fixed point** (`plasticity_2d_vonMises_fixed_point`, id `PLASTICITY_2D_ONECONE_VONMISES_FP`):
sigma <- P_C(sigma − rho (M sigma + q)). Default: 1000 iterations, tol 1e-14.

The step rho is set by `iparam[PLASTICITY_NSN_RHO_STRATEGY]`:

* `..._RHO_STRATEGY_CONSTANT`: rho = `dparam[PLASTICITY_NSN_RHO]` (0 means rho0);
* `..._SPECTRAL_NORM` (and the other spectral variants): rho = rho0;
* `..._ADAPTIVE` (default), as in the friction local iteration solvers:
  * rho is multiplied by 2/3 until rho ||w(sigma+) − w(sigma)|| <= 0.9 ||sigma+ − sigma||, at
    most `iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER]` times (default 10);
  * rho is multiplied by 3/2 when that ratio is below 0.3.

**Semi-smooth Newton** (`plasticity_2d_vonMises_semismooth_newton`, id
`PLASTICITY_2D_ONECONE_VONMISES_NSN`), on Phi(sigma) = sigma − P_C(sigma − rho w(sigma)).
Default: 50 iterations, tol 1e-14.

* Generalized Jacobian: H = I − J_P(y)(I − rho M), with y = sigma − rho w. J_P = I inside C.
  Outside C, J_P is the identity on sigma_N and (R/||y_T||)(I − n n^T) on sigma_T, with
  n = y_T/||y_T||.
* Options:
  * rho strategy: constant, or spectral norm (default);
  * `iparam[PLASTICITY_NSN_LINESEARCH]`: Armijo (default), Goldstein–Price or none, on
    1/2 ||Phi||^2, with `iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER]` steps (default 20);
  * `iparam[PLASTICITY_NSN_HYBRID_STRATEGY]`: any value other than `NO` (the default is one of
    them) falls back to the return mapping when the Newton method does not converge.

### Plane stress: the ball model

**Yield set.** In plane stress (szz = 0), with p = (sxx + syy)/2 and d = (sxx − syy)/2,

    sigma_eq² = sxx² + syy² − sxx syy + 3 sxy² = p² + 3 (d² + sxy²)
              = sigma_N²/2 + (3/2) ||sigma_T||²   <= sigma_y²

This is an ellipsoid with semi-axes sqrt2 sigma_y along N and sqrt(2/3) sigma_y along T. The
hydrostatic part is bounded:

* equibiaxial tension sxx = syy = sigma_y is on the surface (sigma_N = −sqrt2 sigma_y);
* uniaxial tension sxx = sigma_y is on the surface;
* pure shear yields at sxy = sigma_y/sqrt3, as in plane strain.

The flow is normal to the ellipsoid, so w_N ≠ 0: in plane stress the in-plane flow is not
isochoric, because ezz is free.

**Change of variables.** With

    sigma = D y,   D = diag(sqrt2, sqrt(2/3), sqrt(2/3)),

the ellipsoid becomes the ball ||y|| <= sigma_y (D does not depend on sigma_y, so it is the
same for all the elements). The strain rate is mapped with the dual map, w_y = D w, so the
power is preserved. The driver composes D with the Cartesian map (T = P^{-1} D), solves the
problem in y with the hypothesis `PLASTICITY_VONMISES_BALL`, and maps the solution back. The
solution of the problem is the same as with a projection on the ellipsoid in the
sigma-metric: only the metric of the natural map changes (the error is measured in y).

`plasticity_2d_compute_error` on a plane stress problem computes the error in y:
y = D^{-1} sigma, w_y = D w, error of the ball.

**Ball local solvers** (`plasticity_2d_vonmises.c`, selected by the hypothesis of the local
problem, with the same solver ids as the cylinder):

* projection: P_B(y) = y min(1, sigma_y/||y||) (`plasticity_2d_projectionOnVonMisesBall`);
* return mapping (`plasticity_2d_vonMises_ball_return_mapping`):
  * elastic trial y = −M^{-1} q;
  * plastic case: (M + mu I) y = −q with ||y|| = sigma_y, mu = lambda/sigma_y > 0;
  * the scalar equation 1/||y(mu)|| = 1/sigma_y is solved by the same safeguarded Newton as
    for the cylinder (`vonmises_shifted_norm_equation`, here in dimension 3 and without the
    Schur complement on N);
  * this is the classical plane stress closest-point return mapping, in the M-metric;
* fixed point and semi-smooth Newton on the natural map
  (`Plasticity2DVonMisesLocalParameters::ball = 1`): J_P = I inside the ball and
  (sigma_y/||y||)(I − u u^T), u = y/||y||, outside;
* one projection step (`PLASTICITY_2D_ONECONE_ProjectionOnCone`): y <- P_B(y − rho w),
  rho = 1/M_NN.

MixedProblems Von Mises blocks are plane strain only.

## NSGS and local solvers (`plasticity_2d_nsgs.c`, `plasticity_2d_nsgs_generic_adapter.c`)

| internal solver id | Drucker–Prager | Von Mises |
|---|---|---|
| `PLASTICITY_2D_ONECONE_NSN`, `..._NSN_GP`, `..._NSN_GP_HYBRID` (default of the NSGS) | semi-smooth Newton | semi-smooth Newton |
| `PLASTICITY_2D_ONECONE_VONMISES_NSN` | — | semi-smooth Newton |
| `PLASTICITY_2D_ONECONE_ProjectionOnCone` | projection step | one projection step sigma <- P_C(sigma − rho w), rho = 1/M_NN |
| `PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration` | projection with local iterations | fixed point |
| `PLASTICITY_2D_ONECONE_VONMISES_FP` | — | fixed point |
| `PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING` (`PLASTICITY_2D_ONECONE_VONMISES`, deprecated) | — | return mapping |

The Drucker–Prager ids select the Von Mises solver of the same family, so the default options
of the NSGS work for both models. MixedProblems Von Mises blocks use the return mapping.

## Bugs fixed (2026-11-23)

* `plasticity_2d_nsgs_init_wrap` (registry `init` of `PLASTICITY_2D_NSGS`) called
  `plasticity_2d_nsgs_set_default` at each driver call. Whatever the user had set, the
  tolerance was reset to 1e-4 and the internal solver replaced by the default one. This is why
  the NSGS "converged" with an error ~1e-6 for a tolerance of 1e-14 (footing problems). The
  same pattern was fixed in `fc2d_nsgs` (and its parallel variants), the fc3d
  Panagiotopoulos / Tresca / ACLM fixed-point solvers and `fc3d_proximal`.
* `PLASTICITY_2D_NSGS_GENERIC` did not iterate at all: its wrapper passed `info = 0` to
  `nsgs_solve`, which means "trivial case already solved", and
  `test_plasticity_2d_nsgs_generic` did not check the result. The adapter also printed every
  iteration (`verbose = 2` hard coded); it now follows the numerics verbosity.
* `plasticity_2d_compute_error` rejected the Von Mises model.
* (2026-11-25) Python bindings: the problems were freed with `delete` (M, q and the model
  parameters leaked), and M was copied from a C-ordered numpy array into the column-major
  dense storage, i.e. transposed for a non-symmetric M. The problems are now created by
  `plasticityProblem_new` and freed by `plasticityProblem_free`, and M is converted to
  Fortran order. `plasticity_2d_driver` was not exposed, so no Python test solved a
  problem.
* (2026-11-25) The final full-error test of the NSGS (`LIGHT_WITH_FULL_FINAL`, the default)
  compared the error with the incremental tolerance, which the loop may have lowered, instead
  of the required tolerance. The NSGS then returned `info = 1` with a residual below the
  tolerance (seen with the one-step projection local solver on the plane stress problems). The
  friction and rolling-friction NSGS already used `SOLVER_TOL(options)`.
* Trivial-solution tests (the driver and the Newton local solver) returned sigma = 0 whenever
  q_N >= 0; see the correct conditions above.
  * With the default Newton local solver, the NSGS stalled on
    `plasticity_2d_footing_100_theta0.dat` (theta = 0, eta = 0.1), with a cone kept at the
    apex.
  * It now converges with all local solvers (~2300 iterations, error ~2e-12).

## Tests (`numerics/src/Plasticity/test`)

* `plasticity_2d_simple_solve_test`: example1 (theta, theta = 0) and footing_1 with the NSGS
  (default and Newton local solvers); footing_100 (theta, theta = 0.05, theta = 0) with the
  generic NSGS, tolerance 1e-10.
* `plasticity_2d_vonmises_test`:
  * 1000 pseudo-random non-symmetric local problems (condition numbers up to ~1e2):
    * the return mapping has a residual < 2e-16;
    * Newton takes 2.6 iterations on average (all line searches);
    * the fixed point takes 100 (adaptive rho) / 184 (rho = 1/||M||_2) iterations;
    * all agree with the return mapping;
  * NSGS with all local solvers on `plasticity_2d_vonmises_example1.dat` and the footing
    problem;
  * a non-trivial case with q_N >= 0.
* `plasticity_2d_cohesion_test`:
  * M = I analytic cases: tension within the cohesion, apex, sliding;
    sigma = P_{K_eta}(−q + a) − a;
  * a cone with q_N >= 0 but theta < eta (not trivial), with three local solvers;
  * footing problem with c = 1% of the mean normal stress: same solution as the translated
    problem solved by hand.
* `plasticity_2d_cartesian_test` (the map is written again from the definitions):
  * Von Mises with M = I (Voigt, Mandel), pure shear and plane strain uniaxial loading:
    analytic values to ~1e-16;
  * random coupled non-symmetric problems (4 elements, 3 or 4 on the yield surface):
    * models: Drucker–Prager with dilatancy and cohesion, and Von Mises;
    * Mandel and Voigt; dense, sparse-block and sparse M;
    * the Cartesian solution equals the cone solution mapped back (~1e-16);
    * the power is preserved;
    * `plasticity_2d_compute_error` gives ~1e-15;
    * `plasticity_2d_problem_to_cone` gives back the cone problem;
  * footing problem rewritten in Voigt coordinates: same solution (relative difference 5e-11
    at tolerance 1e-12) and same number of iterations.
* `plasticity_2d_plane_stress_test`:
  * ball local solvers on 1000 pseudo-random non-symmetric local problems (959 plastic):
    * the return mapping has a residual < 4e-16;
    * the fixed point (adaptive rho, 22 iterations on average) and the semi-smooth Newton
      (3.4 iterations) agree with it to ~5e-14;
  * manufactured solutions (the plane stress formula and its gradient are written again in
    the test, in Voigt, Mandel and cone coordinates):
    * target stresses: uniaxial tension and compression, equibiaxial tension and compression,
      pure shear, random points of the yield surface, an interior point;
    * a coupled non-symmetric M with a positive definite symmetric part, and
      q = −M s − lambda grad(f²/2)(s), so that the solution is unique and equal to s;
    * dense, sparse-block and sparse M; default, return mapping, fixed point and one-step
      projection local solvers;
    * results: the solution equals s (~1e-16), sigma_eq = sigma_y on the plastic elements,
      and `plasticity_2d_compute_error` gives ~1e-15;
  * `plasticity_2d_problem_change_of_variables`: the ball problem gives s = T y.
* Python, `numerics/python/tests/test_plasticity_2d.py` (44 cases, through
  `plasticity_2d_driver` and `plasticity_2d_compute_error`). The coordinate maps and the plane
  stress formula are written again in the test:
  * creation, invalid sizes, layout of M (C or Fortran order), default options;
  * Drucker–Prager: M = I analytic cases with cohesion, a non-trivial case with q_N >= 0,
    and a coupled problem with dilatancy and cohesion (all the local solvers agree);
  * Von Mises plane strain: pure shear and uniaxial loading (Voigt, Mandel);
  * Von Mises plane stress: the manufactured solutions of the C test (non-symmetric M;
    Voigt, Mandel and cone coordinates; all the Von Mises local solvers), and equibiaxial
    tension (on the plane stress surface, inside the plane strain cylinder);
  * Cartesian Drucker–Prager problems equal to the cone problem mapped back (power
    preserved).
* MixedProblems: `mixed_problem_mohr_coulomb_test` (coupled cones with cohesion, mixed NSGS
  vs Plasticity NSGS).

## Known problems

* `plasticity_2d_footing_10.dat` is in an old format (no model type line) and cannot be read.
* Neither the cohesion, the stress coordinates nor the Von Mises hypothesis are saved in the
  data files.
* Not covered: elasto-plastic plane strain with szz as an unknown; plane stress in
  MixedProblems.
