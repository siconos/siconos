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
| `plasticity_2d_druckerprager.{h,c}` | Drucker–Prager local solvers: projection, projection with local iterations, return mapping, semi-smooth Newton |
| `plasticity_2d_local_newton.{h,c}` | semi-smooth Newton on the natural map of one element (line searches, rho, fallback), shared by Drucker–Prager and Von Mises |
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
* `plasticity_printInFile` saves the coordinates, the Von Mises hypothesis and the cohesion
  (see "Data file format").

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

### Return mapping: exact one-cone solver (`plasticity_2d_DruckerPrager_return_mapping`)

For a 3x3 local matrix M (any, not necessarily symmetric) and any dilatancy theta, the
solver enumerates the three possible cases:

* **elastic**: w = 0, sigma = −M^{-1} q, if it lies in K;
* **apex**: sigma = 0, if q + theta ||q_T|| e_N is in the dual cone;
* **sliding**: sigma on the boundary, sigma = sigma_N (1, eta t) with t = (cos phi, sin phi)
  and sigma_N > 0. The normality rule gives w_T = −lambda t and w_N = lambda (eta − theta),
  lambda >= 0, so

      sigma_N M a(phi) − lambda b(phi) = −q,   a = (1, eta cos phi, eta sin phi),
                                               b = (eta − theta, −cos phi, −sin phi).

  These 3 equations in (sigma_N, lambda, phi) are solvable only if det[M a, b, q] = 0:

  * this determinant is bilinear in a and b, i.e. a trigonometric polynomial of degree 2 in
    phi;
  * with tau = tan(phi/2), (1 + tau²)² det is a quartic in tau, solved by `BIQUADROOTS`
    (`utils/quartic.h`; lower degree when the leading coefficients vanish, phi = pi added);
  * for each real root: least-squares (sigma_N, lambda), then Newton's method on
    (sigma_N, lambda, phi) to round-off.

The solutions are the candidates with a zero natural-map residual (up to round-off).

* For a symmetric positive definite M and theta = eta (associated flow), the solution is
  unique: the closest-point projection of the trial stress −M^{-1} q on K in the M-norm.
* Otherwise (non-associated flow, non-symmetric M) there may be several solutions. The one
  closest to the input stress is returned, for continuity of the Gauss–Seidel iterations.

### Semi-smooth Newton (`plasticity_2d_DruckerPrager_semismooth_newton`)

The method works on the natural map

    Phi(sigma) = sigma − P_K(sigma − rho v),   v = w + theta ||w_T|| e_N,

with the generalized Jacobian H = I − J_P(y)(I − rho J_v):

* J_v = (I + theta e_N [0, n^T]) M, with n = w_T/||w_T|| (0 if w_T = 0);
* J_P is the Jacobian of the projection on K, with the same cases as `soc_projection_3d`:
  * 0 in the polar cone;
  * I inside K;
  * on the boundary, with P = alpha (1, eta u), alpha = (y_N + eta ||y_T||)/(1 + eta²) and
    u = y_T/||y_T||, the derivative of this formula.

The line searches, the rho strategy, the stopping test and the fallback to the return mapping
are those of the Von Mises Newton (`plasticity_2d_local_semismooth_newton`). The existing
`PLASTICITY_2D_ONECONE_NSN*` solvers (Alart–Curnier or natural map, ported from fc3d) remain:
they use a rho per component, static state, no line search by default, and a hybrid
projection loop.

On the random local problems of `plasticity_2d_druckerprager_test`, without fallback, the
Newton method takes 2.0 iterations on average and never fails.

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

## Error and scaling (`plasticity_2d_compute_error`)

The error is scaled as in `fc3d_compute_error`, with the same option indices and values.

For each element e, the quantities are computed in the variables of the solvers (cone
coordinates, cohesion shift, plane stress scaling):

* r_e, the stress;
* v_e, the modified strain rate: w + theta ||w_T|| e_N for Drucker–Prager, w for Von Mises;
* F_e = r_e − P_K(r_e − rho_e v_e), the natural-map residual. rho_e balances the stress and
  the strain rate, which have different units: see "Balancing rho" below.

The error is then:

* scope `PLASTICITY_SCALING_SCOPE_PER_CONE` (default): sqrt(sum_e ||F_e||² / s_e²), where s_e
  is computed from r_e and v_e;
* scope `PLASTICITY_SCALING_SCOPE_PER_VECTOR`: ||F|| / s, where s is computed from the whole
  vectors.

The scale combines the norms of r and v, each safeguarded by max(1, ·):

| option | values |
|---|---|
| `iparam[PLASTICITY_SCALING_NORM]` (24) | `EUCLIDEAN` (default); `NO` (no scaling); `CONE` \|x_N\| + \|\|x_T\|\|; `CONE_ETA` \|x_N\| + eta \|\|x_T\|\| (Von Mises: eta = 1); `CONE_WEIGHTED` sqrt(x_N² + alpha \|\|x_T\|\|²) |
| `iparam[PLASTICITY_SCALING_REDUCTION]` (26) | `MIXED` (default): \|r\| + \|v\|; `MAX`: min(max(1, \|r\|), max(1, \|v\|)), named as in fc3d (its default); `MEAN`: geometric mean |
| `iparam[PLASTICITY_SCALING_SCOPE]` (25) | `PER_CONE` (default), `PER_VECTOR` |
| `iparam[PLASTICITY_SCALING_RHO]` (27) | `ELEMENT` (default): rho_e = 1/\|\|M_ee\|\|_2; `GLOBAL`: rho = 1/max_e \|\|M_ee\|\|_2; `ONE`: rho = 1 (as in fc3d) |
| `dparam[PLASTICITY_SCALING_ALPHA]` (5) | alpha of `CONE_WEIGHTED` |

**Defaults.** `plasticity_2d_nsgs_set_default` (used by `PLASTICITY_2D_NSGS` and
`PLASTICITY_2D_NSGS_GENERIC`) sets `EUCLIDEAN`, `PER_CONE`, `MIXED` and `ELEMENT`, and
`options == NULL` gives the same defaults. The error is then relative to the stress of each
element, and it does not depend on the units of the strain rate. The NSGS passes its options
to the error, so these options are set on the NSGS options. The enum values are those of
fc3d: options created for another solver (zero-initialized) give the fc3d scaling (`MAX`,
rho = 1). `NO` with `ONE` gives the former absolute error. Python: `sn.params.PLASTICITY_SCALING_*`, and
`plasticity_2d_compute_error(problem, stress, w, tol, options=None)`.

**Scaling on the footing problems.** There, |sigma_e| ~ 3e2–5e3 while |w_e| < 4e-4 (M is a
compliance, ~1e-6). The default scaling (min of the safeguarded norms) is then 1, so the
error stays absolute. With the new Drucker–Prager Newton local solver:

| scaling | tol 1e-10 | tol 1e-12 | tol 1e-14 |
|---|---|---|---|
| `NO` (absolute) | ok (error 2e-12) | stalls at 2e-12 | stalls |
| default (`EUCLIDEAN`, `PER_CONE`, `MAX`) | ok (2e-12) | stalls at 2e-12 | stalls |
| `PER_CONE`, `MEAN` | ok (5e-14) | ok, 3016 iterations | stalls at 4e-14 |
| `PER_CONE`, `MIXED` | ok (1e-15) | ok, 3016 iterations | ok, 3636 iterations (8e-16) |
| `PER_VECTOR`, `MIXED` | ok (9e-17) | ok, 3016 iterations | ok, 3636 iterations |

At tolerance 1e-10 the NSGS stops after the same number of iterations (2388) with every
scaling: the stop is decided by the light (incremental, relative) error, and the final
full-error check passes.

**Balancing rho.** The natural map with rho = 1 compares a stress with a strain rate. The
scales of r and v (and the reductions) are then meaningful only if both have comparable
magnitudes, which is not the case for compliance problems (|sigma| ~ 5e3, |w| < 4e-4 on the
footing problems). With rho_e = 1/||M_ee||_2, computed on the diagonal block T^T M_ee T in
the variables of the solvers, rho_e v_e is a stress. The solutions are the same (zero of the
natural map for any rho > 0), and the error no longer depends on the units of the strain
rate: (kappa M, kappa q) gives the same error, which is checked by the test.

The choice of the reduction still matters with rho_e. In each element one of the two
variables usually vanishes: w = 0 for an elastic element (sticking contacts in fc3d). Then
the per-cone `MAX` scale min(max(1, |r_e|), max(1, |rho_e v_e|)) is 1, and the error of the
element stays absolute in stress units (floor eps |sigma_e|). On the footing problem, with
rho_e and the new Drucker–Prager Newton local solver:

| scaling (rho_e) | tol 1e-12 | tol 1e-14 |
|---|---|---|
| `PER_CONE`, `MAX` | stalls at 4e-12 | stalls |
| `PER_CONE`, `MEAN` | ok, 3600 iterations | stalls at 1e-13 |
| `PER_CONE`, `MIXED` | ok, 3323 iterations | 1.06e-14 (floor ~1e-14) |
| `PER_VECTOR`, `MAX` | ok, 3296 iterations | 1.01e-14 |
| `PER_VECTOR`, `MIXED` | ok, 3016 iterations | ok, 3636 iterations (8e-15) |

The objection to `MIXED` (it under-reports the error of the smaller variable when the two
have different magnitudes, `Error_SOC_natural_map.tex`) is removed by rho_e. `PER_CONE` +
`MIXED` + `ELEMENT`, the default since 2026-11-25, gives an error relative to the stress of
each element.

With the default error, all the Drucker–Prager local solvers converge at tolerances 1e-10
and 1e-12 on the footing problems (see the table of the next section). The floor of the
default error is ~1e-14.

**Local errors.** The local solvers measure their error with the same scaling, through
`plasticity_2d_DruckerPrager_local_scaled_error`: ||F|| / max(1, ||r|| + rho ||v||), with
rho = 1/||M_loc||_F. This is used by the old Newton solvers (`plasticity_2d_compute_local_error`)
and by the projection with local iterations. The local tolerance (1e-14 by default) is thus
relative, like the global one. Before, the old Newton solvers used ||F|| / ||q_loc||, which
on compliance problems (|q_loc| ~ 1e-4) could not reach the local tolerance. The projection
with local iterations used the squared absolute residual.

The Frobenius norm is used for rho (local and global errors): it is cheap and within a factor
sqrt(3) of ||M||_2. With a power iteration at every evaluation, the old Newton solvers were
15 times slower.

## Global semi-smooth Newton (`PLASTICITY_2D_NSN`, `plasticity_2d_nsn.c`)

NSGS needs O(cond(M)) sweeps: on the problems of a FEM solid in a quasi-static regime
(`plasticity_2d_vonmises_strip_*.dat`, cond(M) ≳ 1e4) it takes 5e4–7e4 sweeps, and its
stresses are only accurate to about cond(M) × tolerance. `PLASTICITY_2D_NSN` solves the whole
problem by a semi-smooth Newton method. It works on the problem given by
`plasticity_2d_driver` after the change of variables (cone coordinates; Von Mises cylinder or
ball; Drucker–Prager without cohesion).

**Natural map.** For each element e, with w = M σ + q,

    Φ_e(σ) = σ_e − P_{K_e}(σ_e − ρ_e v_e(w)),

v_e = w_e (Von Mises), v_e = w_e + θ_e ‖w_eT‖ e_N (Drucker–Prager, dilatancy). An element of
the generalized Jacobian is

    H = (I − J_P) + J_P R G M,

J_P = blockdiag(J_P,e) the Jacobians of the projections at y_e = σ_e − ρ_e v_e
(`plasticity_2d_vonMises_projection_jacobian`, `plasticity_2d_DruckerPrager_projection_jacobian`,
the same as in the one-element Newton solvers), R = blockdiag(ρ_e I), G = blockdiag(∂v_e/∂w_e).
H has the sparsity of M; the Newton direction H d = −Φ is computed by a sparse LU (`NM_LU_solve`).
H was checked against finite differences of Φ (relative difference 1e-10).

**Globalization** (`iparam[PLASTICITY_NSN_GLOBALIZATION]`):

* `WATCHDOG` (default): full Newton steps σ ← σ + d. The active set (plastic elements) is
  corrected for all the elements at once, like a primal-dual active set method; the residual
  is not monotone. After `PLASTICITY_NSN_WATCHDOG_ITERATIONS` (50) full steps without
  convergence, a non-finite iterate or a singular system, the iteration restarts from the best
  iterate with the merit function below.
* `MERIT`: merit function only, chosen automatically:
  * *QP objective* when the problem is a convex QP — associated flow (Von Mises, Drucker–Prager
    with θ = 0) and M symmetric (`‖M − Mᵀ‖∞ ≤ 1e-12 ‖M‖∞`): f(σ) = ½ σᵀMσ + qᵀσ over K, feasible
    iterates (the initial guess is projected), Armijo search along the projected arc
    σ(t) = P_K(σ + t d) (f(σ(t)) ≤ f(σ) + 1e-4 wᵀ(σ(t) − σ) + 16 ε |f(σ)|, the last term
    accepting the decreases below the round-off of f), projected gradient arc P_K(σ − t R w)
    when the Newton arc gives no decrease (projected Newton, Bertsekas);
  * *natural map* ½‖Φ‖² otherwise (dilatancy, non-symmetric M): Armijo line search,
    Levenberg–Marquardt direction `(HᵀH + λI) d = −HᵀΦ`, λ = c‖Φ‖, when the iteration stalls or
    H is singular (`nonsmooth_newton_levenberg_marquardt.h`).
* `NONE`: full steps only (returns the best iterate).

**Choice of ρ** (`iparam[PLASTICITY_NSN_RHO_STRATEGY]`, `dparam[PLASTICITY_NSN_RHO_FACTOR]`).
ρ_e decides which elements a full Newton step treats as plastic (through y_e), so it shapes
the path of the full steps. Default: ρ_e = c ‖M_ee⁻¹‖₂ = c/λ_min(M_ee) with c = 10 for the
full steps and c = 1 in the merit phase (`INVERSE_SPECTRAL_NORM`; `SPECTRAL_NORM`:
ρ_e = c/‖M_ee‖₂). Experiments on the FEM strip problems (60 to 1446 elements) and the footing
problems, number of full steps from a zero initial guess:

| ρ_e | strip 1446, step 10 | strip 636, step 18 | strip 60, step 100 | footing 100 |
|---|---|---|---|---|
| 1/‖M_ee‖₂ | – | 73 | 19 | 11 |
| ‖M_ee⁻¹‖₂ | 80 | 19 | 11 | 12 |
| 10 ‖M_ee⁻¹‖₂ (default) | 23 | 12 | 9 | 9 |
| 30 ‖M_ee⁻¹‖₂ | 18 | 11 | 7 | 12 |
| 100 ‖M_ee⁻¹‖₂ | 15 | 11 | 7 | diverges |
| ‖(M⁻¹)_ee‖₂ | 8 | 6 | 4 | diverges |

The full steps converge fastest for ρ in a window (too large: no convergence). The diagonal
block of M⁻¹ (the effective compliance of the element in the coupled problem, 5e3–2e4 times
‖M_ee⁻¹‖ on the strip) is at the edge of this window and fails on several problems; c = 10 is
robust on all the problems tried. With c = 10 the one-cone `footing_1` problem cycles and is
solved by the watchdog (50 full steps, then 2 merit iterations); the merit phase converges
faster with c = 1 (strip step 10: 33 iterations, 196 with c = 10).

**Options.** `iparam[SICONOS_IPARAM_MAX_ITER]` (100), `dparam[SICONOS_DPARAM_TOL]` (1e-10,
`plasticity_2d_compute_error`, the error of the NSGS solvers),
`iparam[PLASTICITY_NSN_GLOBALIZATION]` (`WATCHDOG`), `iparam[PLASTICITY_NSN_RHO_STRATEGY]`
(`INVERSE_SPECTRAL_NORM`), `dparam[PLASTICITY_NSN_RHO_FACTOR]` (10),
`iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER]` (30), `iparam[PLASTICITY_NSN_LM_FALLBACK]` (1),
`dparam[PLASTICITY_NSN_LM_COEF]` (1).

**Results** (macOS arm64, zero initial guess; NSGS tolerance 1e-8, NSN 1e-10, default options):

| problem | NSGS sweeps | NSN iterations |
|---|---|---|
| strip step 1 (elastic) | 72 894 | 1 |
| strip step 4 (1 plastic element) | 72 365 | 4 |
| strip step 10 (18 plastic) | 53 216 | 7 |
| strip step 100 (23 plastic) | 52 575 | 9 |
| footing 100, θ > 0 (natural map) | | 9 |
| footing 100, θ = 0 | | 9 |
| Von Mises footing 100 | | 5 |

NSN and NSGS stresses differ by 4e-5 relative on the strip: the error of NSGS
(cond(M) × 1e-8). In the FEM simulation of the strip (100 steps, warm start from the previous
step, tolerance 1e-10; meshes of pytripopfem):

| mesh (T3 elements, stresses) | NSN iterations per plastic step (mean / max) | time | final reaction |
|---|---|---|---|
| 60, 180 | 8 / 9 | 0.33 s | 56.51 |
| 636, 1 908 | 13 / 18 | 24 s | 53.33 |
| 1 446, 4 338 | 19 / 30 | 167 s | 52.94 |

(NSGS: about 100 s on the 60-element mesh.) The number of iterations grows slowly with the
mesh. On the 1446-element mesh, a Newton iteration costs ~54 ms (new sparse LU of H at each
iteration): ~100 s of the 167 s; the rest is the kernel (dense iteration matrix and its LU at
the initialization, assembly of W and H at each step, see `FEM_PLASTICITY_SCHEME.md`).

**Behaviour and limits.**

* With the merit function alone, the iteration first identifies the plastic elements with
  short projected-arc steps (the Newton direction is computed for a wrong active set and
  overshoots), then converges fast; this first phase grew with the mesh (636 elements: up to
  470 iterations), hence the full steps by default.
* The error stagnates at ~1e-13 on the strip problems (round-off amplified by cond(M)):
  tolerances below ~1e-12 are not reached reliably; 1e-10 is the practical value.
* Possible improvements: reuse the symbolic factorization of H (same pattern at each
  iteration), a better ρ (between ‖M_ee⁻¹‖ and ‖(M⁻¹)_ee‖), a theory of the window of ρ.

## NSGS and local solvers (`plasticity_2d_nsgs.c`, `plasticity_2d_nsgs_generic_adapter.c`)

| internal solver id | Drucker–Prager | Von Mises |
|---|---|---|
| `PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN` (default of the NSGS since 2026-11-25) | semi-smooth Newton (shared implementation) | semi-smooth Newton |
| `PLASTICITY_2D_ONECONE_NSN`, `..._NSN_GP`, `..._NSN_GP_HYBRID` (`NSN_GP_HYBRID` was the default) | semi-smooth Newton (Alart–Curnier / natural map, fc3d port) | semi-smooth Newton |
| `PLASTICITY_2D_ONECONE_VONMISES_NSN` | semi-smooth Newton (shared implementation) | semi-smooth Newton |
| `PLASTICITY_2D_ONECONE_ProjectionOnCone` | projection step | one projection step sigma <- P_C(sigma − rho w), rho = 1/M_NN |
| `PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration` | projection with local iterations | fixed point |
| `PLASTICITY_2D_ONECONE_VONMISES_FP` | projection with local iterations | fixed point |
| `PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING` | return mapping | return mapping |
| `PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING` (`PLASTICITY_2D_ONECONE_VONMISES`, deprecated) | return mapping | return mapping |

Each id selects the solver of the same family for the other model, so the default options of
the NSGS work for both models.

On the footing problems (100 cones, default error), all the Drucker–Prager local solvers
converge in the same number of NSGS iterations:

| local solver | tol 1e-10 | tol 1e-12 |
|---|---|---|
| `DRUCKERPRAGER_NSN` (default) | 0.68 s | 0.83 s |
| `DRUCKERPRAGER_RETURN_MAPPING` | 0.80 s | 0.98 s |
| `ONECONE_NSN` | 0.75 s | 0.93 s |
| `ProjectionOnConeWithLocalIteration` | 0.73 s | 0.79 s |
| `ONECONE_NSN_GP`, `ONECONE_NSN_GP_HYBRID` | 4.0 s | 5.0 s |
| `ProjectionOnCone` (one step) | 0.08 s | |

Iterations on footing_100: 2655 at 1e-10, 3275 at 1e-12; footing_100_theta0 is similar.
Before the scaled local error, `NSN_GP_HYBRID` took 6.5–11 s and stalled on
footing_100_theta0 at 1e-12.

MixedProblems Von Mises blocks use the return mapping.

With the former absolute error (or the fc3d scaling, see "Error and scaling") and
|sigma| ~ 5e3, a tolerance of 1e-12 is at round-off (relative 5e-16): the NSGS stalls at
~2e-12 with every local solver.

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
* (2026-11-25) The error of `plasticity_2d_compute_error` was absolute, whatever the scale of
  the problem (the `norm` argument was unused). It is now scaled as in fc3d (see "Error and
  scaling"); `PLASTICITY_SCALING_NO` gives the former error.
* (2026-11-25) Generic NSGS (`NonSmoothSolvers/nsgs_generic.h`, `nsgs_generic_instrumented.h`,
  shared with the friction generic solvers): after the final check of the full error against
  the user tolerance, info was set again by comparing the error with the working tolerance,
  which may have been lowered during the iterations. The solver then returned `info = 1`
  with an error below the user tolerance. With a full error, the final check now decides.
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

## Data file format (`plasticity_printInFile`, `plasticity_newFromFile`)

    dimension
    numberOfCones
    model_type                       (PlasticityModelType)
    M                                (NM_write_in_file)
    q                                (dimension * numberOfCones values)
    eta, then theta (Drucker-Prager) | sigma_y (Von Mises)   (one line each)
    [stress_coordinates <int>]       (if not PLASTICITY_STRESS_COORDINATES_CONE)
    [vonmises_hypothesis <int>]      (if not PLASTICITY_VONMISES_PLANE_STRAIN)
    [cohesion <numberOfCones values>] (Drucker-Prager with cohesion)

The keyword lines (November 2026) are optional and written only when different from the
defaults: the files of default problems are unchanged, and older files are read as before
(defaults). Before, these three fields were lost (a plane stress / Voigt problem was read
back as a plane strain problem in cone coordinates). Sparse matrices are written with 17
significant digits (`CSparseMatrix_print_in_file`, `%.17g`, was `%g`: 6 digits).
`plasticity_2d_read_write_test` checks the round trip of a plane stress Voigt Von Mises
problem and of a Drucker–Prager problem with cohesion in Mandel coordinates.

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
* `plasticity_2d_druckerprager_test`:
  * return mapping on 3000 manufactured solutions (non-symmetric and symmetric M,
    theta in {0, eta/2, eta, 1.5 eta}):
    * elastic: sigma inside K, w = 0;
    * apex: sigma = 0, v strictly in the dual cone;
    * sliding: sigma on K, w from the flow rule;
    * all recovered to ~1e-15;
  * 2000 random local problems:
    * the return mapping has a residual < 4e-15;
    * the Newton method without fallback always converges (2.0 iterations);
    * both agree where the solution is unique (theta = eta, M symmetric);
  * NSGS on the three footing problems with five local solvers: same solution (see the
    timings above).
* `plasticity_2d_compute_error_test`:
  * the unscaled, default, MIXED per-cone and MIXED per-vector errors, computed again in the
    test (Drucker–Prager with dilatancy and cohesion, Von Mises);
  * homogeneity: with (M, 1e4 q) and 1e4 sigma, the absolute error is multiplied by 1e4 and
    the MIXED error is unchanged;
  * rho_e: the error is unchanged for (kappa M, kappa q), kappa = 1e-6 (per element and
    global rho); with rho = 1 it changes;
  * footing problem with rho_e and the MIXED scaling: converges at tolerance 1e-13, with an
    absolute error at round-off (2.5e-12).
* Python, `numerics/python/tests/test_plasticity_2d.py` (49 cases, including the error
  scaling; through
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
* `plasticity_2d_nsn_test`: `PLASTICITY_2D_NSN` (tolerance 1e-10) on the Drucker–Prager
  footing problems (θ > 0 and θ = 0), the Von Mises footing, the small examples and the FEM
  strip (both globalizations), and on the FEM footing problems
  `plasticity_2d_druckerprager_footing_fem_40.dat` (associated) and
  `plasticity_2d_druckerprager_footing_fem_psi10_40.dat` (dilatancy angle 10°, θ = 0.249):
  siconos-tutorials `mechanics/FEM/Plasticity/Footing`, 639 T3, plane strain Mohr–Coulomb
  (η = sin 25°, cohesion √2 c cos φ, Voigt coordinates), step 40, written by `setProblemDump`
  (before the FEM stresses switched to the Mandel notation, 2026-11-26: the problems dumped
  now declare `stress_coordinates` Mandel; these files keep the Voigt ones).
  NSN: 12 and 13 iterations (default globalization; the merit function alone does not
  converge in 300 iterations; NSGS: error 3e-4 after 2e4 sweeps). Python: NSN = NSGS on a coupled non-symmetric Drucker–Prager problem with
  dilatancy and cohesion, NSN = exact solution of the manufactured plane stress problem.
* `plasticity_2d_strip_test`: NSGS, generic NSGS (tolerance 1e-8) and NSN (1e-10) on the problems of the
  FEM strip with a circular hole (siconos-tutorials `mechanics/FEM/Plasticity/StripWithHole`,
  60 T3 elements, plane stress, Voigt coordinates (dumped before the switch of the FEM to the
  Mandel notation), yield stress 210),
  `data/plasticity_2d_vonmises_strip_<n>.dat` for the time steps n = 1 (elastic), 4 (first
  plastic element), 10 (18 plastic elements) and 100 (plateau, 23 plastic elements), written
  by `fem::nonsmooth_formulations::Plasticity::setProblemDump`. Checks: convergence, error
  recomputed by `plasticity_2d_compute_error`, Von Mises stress ≤ yield stress, same stresses
  for the two NSGS, NSN = NSGS up to 1e-3 (accuracy of NSGS). These problems come from a quasi-static regime where M is
  ill-conditioned (condition number ≳ 1e4): NSGS needs 1e4–1e5 sweeps; they are the
  reference problems for global solvers.
* MixedProblems: `mixed_problem_mohr_coulomb_test` (coupled cones with cohesion, mixed NSGS
  vs Plasticity NSGS).

## Known problems

* `plasticity_2d_footing_10.dat` is in an old format (no model type line) and cannot be read.
* `PLASTICITY_2D_ONECONE_NSN_GP` (damped Newton, Goldstein–Price line search, rho per
  component) and the hybrid solver built on it: on the footing problems, the scaled local
  error does not go below ~1e-12 on some cones, so the local solves run to their iteration
  limit with local tolerances <= 1e-12. They take 4–5 s, against 0.44 s with a local
  tolerance of 1e-10. `ONECONE_NSN` (no line search), with the same functions, reaches 1e-14.
  Not investigated further: the default local solver is now `DRUCKERPRAGER_NSN`.
* Not covered: elasto-plastic plane strain with szz as an unknown; plane stress in
  MixedProblems.
