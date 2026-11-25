# Energy, work of the forces, FEM node velocities and FEM mass (vnative)

What the mechanics runner writes, on the vnative backend, for the energy
balance of a simulation and for the linear FEM bodies: the `energy_work`
table (run option `output_energy_work`), the node velocities
`fem_velocities_<ds_id>`, the mass matrix `fem_mass_<ds_id>` and the dof map
`fem_dof_indices_<ds_id>`. Read beside
[paraview-vtkhdf-export.md](paraview-vtkhdf-export.md) for the other outputs.

Code: `vkernel/src/siconos/vkernel/io/io.hpp` (`kinetic_energies`,
`work_of_forces`, `fem_velocities`, `mass_matrices`, `contact_work`),
`vkernel/src/siconos/vkernel/py/nonos/bridge.py` (`MechanicsIO`),
`io/python/mechanics_run.py` (`output_energy_and_work`,
`output_fem_velocities`, `output_fem_mass`), `io/python/format_v1.json`.
Tests: `vkernel/src/siconos/vkernel/tests/test_energy_work.py`
(`ctest -R vkernel_python_energy_work`).

## 1. The datasets

| dataset (under `/data`)   | when                                  | layout                                                                                               | dtype |
|---------------------------|---------------------------------------|------------------------------------------------------------------------------------------------------|-------|
| `energy_work`             | each output step, `output_energy_work` | `t, kinetic, force_work, w_n, w_t, w_n_theta, w_t_theta, w_n_neg, w_t_neg`                            | f4    |
| `cf_work`                 | each output step, `output_contact_work` | `t, interaction_id, w_n, w_t, w_n_theta, w_t_theta, status, w_n_pos` (one row per active contact)   | f4    |
| `fem_velocities_<ds_id>`  | each output step, always (vnative FEM) | `t, ds_id`, then the dofs in spatial order, as `fem_displacements_<ds_id>` (attrs `dof_count`, `spatial_order`) | f8    |
| `fem_displacements_<ds_id>` | each output step, always            | `t, ds_id, u` in spatial order                                                                       | f4 (the format says f8: see 5) |
| `fem_mass_<ds_id>`        | once, at initialisation               | group, CSC: `outer`, `inner`, `values`; attrs `format` = "csc", `rows`, `cols`, `nnz`, `mass` = "consistent" | i8, i8, f8 |
| `fem_stiffness_<ds_id>`   | once, at initialisation               | the same CSC group, without `mass`                                                                   | i8, i8, f8 |
| `fem_dof_indices_<ds_id>` | once, at initialisation               | `dofs[i]`: the dof of M and K of column `2 + i` of the spatial-order tables                           | i8    |

- The FEM tables have no row for the initial state: their first row is the
  first step. The kinetic energy and the works likewise.
- M is the consistent mass of the FEM model (Gauss quadrature of N N^T),
  the matrix the integrator uses; it is not lumped. It is constant: the FEM
  is linear, nothing changes `mass_matrix` after the import. A Dirichlet
  condition does not modify M or K: it blocks velocity components.
- M and K are in the dof numbering of the FEM model; u and v in spatial
  order. `fem_dof_indices_<ds_id>` maps one to the other (see 4).
- On a cold restart (running again to a later T), the per-step tables
  append and the once-only groups are kept (their metadata checked).

## 2. The columns of `energy_work`

For the step from t(k) to t(k+1) = t, with v, q at both levels, theta and h
of the vkernel integrator (the `theta_dem` of an `ExplicitDEMOSI`, not the
`theta` run option):

- `kinetic` = sum over the systems of 1/2 v(k+1)^T M v(k+1), in J. A rigid
  body: M = diag(m, m, I), the rotation included. A FEM body: the sparse M.
- `force_work` = sum of h f_theta . v_theta, in J, x_theta = theta x(k+1) +
  (1 - theta) x(k), f the *total* forces, as `MoreauJeanOSI::computeWorkForces`
  of the AOS kernel: fext for a rigid body (gravity), fext - K q for a FEM
  body. It therefore contains the elastic work of the FEM bodies: over a
  step, at theta = 1/2, its FEM part is f^T dq - d(1/2 q^T K q).
- `w_n`, `w_t`: the works of the contact *impulses* over the step, J,
  lambda . (ydot(k) + ydot(k+1)) / 2, normal and tangential, summed over
  the active contacts (rigid and FEM contacts). `w_n_theta`, `w_t_theta`:
  the same with the theta weight. `w_n_neg`, `w_t_neg`: the sums of the
  negative ones.

These are works over *one* step, written at the output steps only: their
cumulative sums (`siconos_energy`) are the works since the start only with
`output_frequency` 1.

## 3. The balance they close

The step of the vkernel is M (v(k+1) - v(k)) = h f_theta + H^T lambda (the
FEM free flight included: A = M + h^2 theta^2 K is that equation solved
for v(k+1)). With v_1/2 = (v(k) + v(k+1)) / 2:

    K(k+1) - K(k) = h f_theta . v_1/2 + lambda . H v_1/2
                  = force_work - (theta - 1/2) h f_theta . (v(k+1) - v(k)) + w_n + w_t

- `force_work` equals the change of K only at theta = 1/2. The difference,
  (theta - 1/2) h f_theta . dv, is the numerical dissipation of the
  theta-method, not a rounding: for a disk in free fall, (theta - 1/2) m g^2
  h^2 per step (theta = 1: 1/2 m g^2 h^2).
- `w_n + w_t` is lambda . H v_1/2 when H is the same at k and k+1 (ydot(k)
  is computed with the H of step k): exact for a disk on a segment, O(h) for
  a turning normal.
- Measured (`test_energy_work.py`): free fall, theta 1 and 0.50001, the
  formula above to 2e-6 J; a disk dropped on a segment, e = 0 and 0.5,
  theta 1/2 and 0.50001, the per-step balance to 1e-5 of K; a free FEM
  square under a constant nodal force, K + 1/2 u^T K u - f^T u constant to
  4e-5 J (scale 1.4e3 J, the float32 of the displacements) at theta = 1/2,
  decaying at theta = 0.6.
- The same rigid scene on the AOS backend (bullet): `kinetic` and
  `force_work` identical in free flight; the contact times differ (see 5).

## 4. Reading them

```python
import h5py
import numpy as np
import scipy.sparse as sp


def csc(group):
    shape = (int(group.attrs["rows"]), int(group.attrs["cols"]))
    return sp.csc_array(
        (group["values"][:], group["inner"][:], group["outer"][:]), shape
    )


with h5py.File("bench_fem_disks.hdf5", "r") as f:
    d = f["data"]
    ew = d["energy_work"][:]
    t, kinetic, force_work, w_n, w_t = ew[:, 0], ew[:, 1], ew[:, 2], ew[:, 3], ew[:, 4]
    for name in d:
        if not name.startswith("fem_mass_"):
            continue
        ds_id = name.removeprefix("fem_mass_")
        M = csc(d[name])
        K = csc(d[f"fem_stiffness_{ds_id}"])
        dofs = d[f"fem_dof_indices_{ds_id}"][:]
        rows_v = d[f"fem_velocities_{ds_id}"][:]
        rows_u = d[f"fem_displacements_{ds_id}"][:]
        v = np.zeros((len(rows_v), M.shape[0]))
        u = np.zeros((len(rows_u), K.shape[0]))
        v[:, dofs] = rows_v[:, 2:]  # spatial order -> dofs of M
        u[:, dofs] = rows_u[:, 2:]
        kinetic_fem = 0.5 * np.einsum("ij,ij->i", v, (M @ v.T).T)
        elastic = 0.5 * np.einsum("ij,ij->i", u, (K @ u.T).T)
        print(ds_id, M.shape, kinetic_fem.max(), elastic.max())
```

## 5. Limits and known issues

- Clamped FEM bodies: the energy_work columns first showed that the step
  was wrong. A Dirichlet condition zeroed the blocked velocities *after* the
  unconstrained solve of the FEM free flight, instead of a solve with the
  operator of the free dofs: the free dofs next to the clamp got spurious
  forces, the clamped square drifted by 47 J at theta = 1/2, and the elastic
  energy of `bench_fem_disks.py` reached 587 J with an input of about 1 J.
  Fixed: the FEM operators are factored and assembled with their Dirichlet
  conditions (`algebra::constrain_dofs`); the clamped square now keeps its
  energy to 3e-6 J, `bench_fem_disks.py` closes its balance to 3e-6 under MJ
  and GOSI, and `fem_t3` matches the AOS kernel's `T3_square_200` to 1.3e-14
  over the whole run (the divergence after the impact, taken for chaos, was
  this defect).
- A cold restart restarts a FEM body from its last output state
  (`fem_displacements_`, `fem_velocities_`, `fem_dof_indices_`); the state
  is read from the float32 displacements: the restarted trajectory is the
  one of a single run to 9e-8 (velocities, scale 2). A file without
  `fem_velocities_` (written before) restarts the FEM bodies from their
  initial state, with a warning.
- `energy_work` and `cf_work` are float32 (format_v1): the balance read
  from the file is good to ~1e-7 of K.
- `fem_displacements_`, `fem_epsilon_`, `fem_sigma_` are declared f8 in
  format_v1.json and written f4 (the h5py default); `fem_velocities_` is
  written f8, as declared.
- AOS backend, bullet, 2D: a disk of radius 1 on a Segment rests at y =
  1.636 (contact point 0.636 above the segment); radius 0.1 at y = 0.3 is in
  contact from the first step. Not investigated.

## 6. Cost

Measured on `bench_fem_disks.py` (`mesh_N5`, 72 dofs, 50 disks, CPU,
step 0.50 ms + 2.25 ms of `updateInteractions`), per output step: the C++
calls 5 to 8 us each; `output_energy_and_work` 0.35 ms;
`output_fem_velocities` 1.13 ms (`output_displacements` 0.84 ms), most of it
the h5py row write; `output_fem_mass` 0.9 ms once. File: 457 B per output
step for `fem_velocities_` (gzip), 13.8 KB once for `fem_mass_`.

10^4 disks dropped on a segment, CUDA build (cuNSearch, `OMP_NUM_THREADS=1`),
20 steps, two runs each: `kineticEnergies` 0.035 ms, `workOfForces` 0.045
ms, `output_energy_and_work` 0.71 ms per output step; `output_results` 12.3
ms with `output_energy_work` against 11.6 ms without, for a step of 5.3 ms
(`computeOneStep`) + 4.8 ms (`updateInteractions`, median). A per-body
Python loop over the same 10^4 bodies (two pybind calls each) costs 32 ms.
The energy output is off by default: no row, no cost (the empty
`energy_work` table is created in every file, as before).
