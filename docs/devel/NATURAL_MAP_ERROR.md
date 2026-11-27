# Projections on convex sets and the natural-map error

Code: `numerics/src/tools/convex_sets/`. Theory of the scaled error on second-order cones:
`docs/devel/Error_SOC_natural_map.tex`.

## 1. Why

Most problems of numerics are complementarity problems over a product of closed convex sets
`K = K_0 x ... x K_{m-1}`:

| problem                         | set of one block                         | dual variable v (before rho)            |
|---------------------------------|------------------------------------------|-----------------------------------------|
| LCP                             | orthant R+                               | w = Mz + q                              |
| Relay                           | box [lb, ub]                             | w                                       |
| FC2D / FC3D, GFC2D / GFC3D      | friction cone K_mu (2D, 3D)              | u + mu ‖u_T‖ e_N (De Saxcé)             |
| Tresca                          | R+ x B(0, R)                             | u                                       |
| Rolling friction (2D, 3D)       | rolling cone (dim 3, 5)                  | u + (mu ‖u_T‖ + mu_r ‖u_R‖) e_N         |
| Plasticity Drucker-Prager       | cone K_eta (cone coordinates)            | w + theta ‖w_T‖ e_N (dilatancy)         |
| Plasticity Von Mises (plane strain) | cylinder R x B(0, sqrt(2/3) sigma_y) | w                                       |
| Plasticity Von Mises (plane stress) | ball B(0, sigma_y) (scaled variables) | D w                                    |
| MixedProblem                    | product of the above                     | per block                               |

Their errors are all built on the natural map `F(r, v) = r - P_K(r - v)` and a scaling. Before
this refactoring each problem had its own copies of the projection (2D cone written twice in
fc2d, 3D cone inlined in fc3d, Von Mises cylinder and ball in Plasticity, box in Relay), of the
norm/reduction dispatch (fc2d, fc3d, plasticity: three identical 4 x 3 switches, plus two
per-cone variants) and of the per-cone / per-vector loops. They now differ only by the set `K`
and by the computation of `v`, which stays in each problem.

## 2. Layout

| file                         | content                                                               |
|------------------------------|-----------------------------------------------------------------------|
| `projection_orthant.{h,c}`   | `orthant_projection`, `box_projection`, diagonal Jacobians            |
| `projection_ball.{h,c}`      | `ball_projection`, `cylinder_projection` (R x B), `positive_cylinder_projection` (R+ x B), Jacobians |
| `projection_soc.{h,c}`       | `soc_projection_2d/_3d`, `soc_projection` (any dim), `soc_dual_projection_3d`, `soc_projection_jacobian_3d`, `soc_projection_3d_stack` |
| `projection_rolling_cone.{h,c}` | `rolling_cone_projection` (5D), `rolling_cone_2d_projection`, `rolling_cone_projection_jacobian` |
| `convex_set.{h,c}`           | `ConvexSet` (type, block dim, number of blocks, per-block parameters p1, p2), `convex_set_project_block` (inline dispatch), `convex_set_projection_jacobian_block`, `convex_set_natural_map(_block)` |
| `natural_map_error.{h,c}`    | scaling options (norms, scopes, reductions), `natural_map_error_norm/_scale`, `natural_map_error` and the accumulate/finalize API for products of sets |

Renamed (no compatibility headers; nothing outside numerics used them):

| old (`linalg/`)                           | new (`tools/convex_sets/`)                    |
|-------------------------------------------|-----------------------------------------------|
| `projectionOnCone(r, mu)`                 | `soc_projection_3d(r, mu)`                    |
| `projectionOnDualCone(u, mu)`             | `soc_dual_projection_3d(u, mu)`               |
| `projectionOnSecondOrderCone(r, mu, n)`   | `soc_projection(r, mu, n)` (now any n; was n = 3 only, exit otherwise) |
| `subdifferentialProjectionOnCone(H, r, mu)` | `soc_projection_jacobian_3d(H, r, mu)`      |
| `projectionOnConeVector(nc, r, mu, work)` | `soc_projection_3d_stack(nc, r, mu, work)`    |
| `naturalMapOnCone`, `naturalMapOnConeVector` | `convex_set_natural_map_block`, `convex_set_natural_map` (v = rho u given by the caller) |
| `projectionOnDisk(r, R)`                  | `ball_projection(r, 2, R)`                    |
| `projectionOnCylinder(r, R)`              | `positive_cylinder_projection(r, 3, R)`       |
| `projectionOnGeneralCylinder(r, R, n)` (empty stub) | `positive_cylinder_projection(r, n, R)` |
| `projectionOnRollingCone`, `projectionOn2DRollingCone`, `projectionOnDualRollingCone`, `subdifferentialProjectionOnRollingCone`, `display_status_rolling_cone` | `rolling_cone_projection`, `rolling_cone_2d_projection`, `rolling_cone_dual_projection`, `rolling_cone_projection_jacobian`, `rolling_cone_projection_status_display` |
| `PROJCONE_DUAL/INSIDE/BOUNDARY`           | `SOC_PROJECTION_POLAR/INSIDE/BOUNDARY` (same values) |
| `PROJRCONE_*`                             | `ROLLING_CONE_PROJECTION_*`                   |
| `SOC_scaling.h`: `SOC_norm_*`, `SOC_scaling_*`, `SOC_geom_scaling_*`, `SOC_mixed_scaling_*` (16 functions) | `natural_map_error_norm(norm, x, m, dim, mu, alpha)`, `natural_map_error_scale(norm, reduction, x, y, m, dim, mu, alpha)` |
| `project_on_box` (Relay_Solvers.h)        | `box_projection(x, n, lb, ub)`                |
| `plasticity_2d_scaling`                   | `natural_map_error_scale`                     |

The Plasticity functions `plasticity_2d_projectionOnVonMises(Ball)`,
`plasticity_2d_vonMises_projection_jacobian` and
`plasticity_2d_DruckerPrager_projection_jacobian` are kept as thin wrappers (Von Mises radius,
return values 1 inside / 0 outside).

## 3. The error

For a stack `K` of `m` blocks, `r` and `v` (v already modified and multiplied by rho):

```
F_i = r_i - P_{K_i}(r_i - v_i)
NO           : ‖F‖
PER_BLOCK    : sqrt( sum_i ‖F_i‖^2 / sigma(r_i, v_i)^2 )
PER_VECTOR   : ‖F‖ / sigma(r, v)
sigma(x, y)  : MAX   min(max(1, ‖x‖), max(1, ‖y‖))     = 1 / max of the two relative errors
               MEAN  sqrt(max(1, ‖x‖) max(1, ‖y‖))
               MIXED max(1, ‖x‖ + ‖y‖)
block norms  : EUCLIDEAN ‖x‖, CONE |x_N| + ‖x_T‖, CONE_MU |x_N| + mu ‖x_T‖ (mu = p1 for
               cones, 1 for the other sets), CONE_WEIGHTED sqrt(x_N^2 + alpha ‖x_T‖^2)
stack norms  : sum of block norms (CONE, CONE_MU), sqrt of sum of squares (EUCLIDEAN, WEIGHTED)
```

Options: the same indices for all the problems, `iparam[24]` norm, `iparam[25]` scope,
`iparam[26]` reduction, `dparam[5]` alpha (`NATURAL_MAP_ERROR_*`; the
`SICONOS_FRICTION_3D_SCALING_*`, `SICONOS_FRICTION_2D_SCALING_*` and `PLASTICITY_SCALING_*`
enumerators are now aliases with unchanged values). `natural_map_error_scaling(options,
defaults)` reads them.

A product of sets of different types (MixedProblem, MLCP) uses one `NaturalMapErrorSum`:

```c
NaturalMapErrorSum sum;
natural_map_error_init(&sum);
natural_map_error_add(&scaling, &K_friction, r_f, v_f, &sum);
natural_map_error_add(&scaling, &K_plasticity, r_p, v_p, &sum);
double e = natural_map_error_value(&scaling, &sum);
```

The norms are accumulated (sums or sums of squares) so that PER_VECTOR is the scale of the
whole product.

## 4. Behaviour changes

The refactoring keeps the error definitions. Differences:

- round-off: `r - (u + mu ‖u_T‖ e_N)` instead of `(r - u) - mu ‖u_T‖`; stack norms by
  accumulated squares instead of `cblas_dnrm2`; Von Mises `sqrt` instead of `hypot`;
- `positive_cylinder_projection` (Tresca, former `projectionOnCylinder`): for `r_N < 0` the old
  code put `r_T` on the boundary even when `‖r_T‖ < R`; the projection on `R+ x B` leaves it
  unchanged (bug fix);
- `projectionOnGeneralCylinder` (SOCLCP projection solver) was an empty function: it now
  projects;
- `projectionOnSecondOrderCone` exited for a dimension other than 3: `soc_projection` handles
  any dimension (SOCLCP);
- `natural_map_error_scale` with the MIXED reduction is safeguarded (`max(1, .)`); the former
  `SOC_mixed_scaling_*` returned the raw sum (the callers safeguarded it).

Not changed: the MLCP error (`mlcp_compute_error`) is not a natural map (sum of `z_i w_i` and
negative parts relative to `1 + ‖q‖`); the global (GFC3D, global rolling) primal/dual residuals
and their `max(‖Mv‖, ‖Hr‖, ‖q‖)` scalings; the relative normalizations by `‖q‖` done by
`lcp_compute_error`, `relay_compute_error`, `rolling_friction_*_compute_error`,
`mixedProblem_compute_error`.

## 5. Adding a set

1. Write the projection (and its Jacobian) in a `projection_<set>.{h,c}` file.
2. Add a `CONVEX_SET_<SET>` type and its case in `convex_set_project_block`,
   `convex_set_projection_jacobian_block`, and in `cone_coefficient` of `natural_map_error.c`
   if the CONE_MU norm should weight it.
3. In the problem: build `v`, a `ConvexSet`, call `natural_map_error`.
