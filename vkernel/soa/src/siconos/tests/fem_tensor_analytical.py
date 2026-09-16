"""
Analytical regression check for FEM strain/stress tensor computation.

This script imposes a synthetic, globally affine displacement field
u(x, y) = (a*x + b*y, c*x + d*y) on the nodes of the FEM mesh and verifies
that every T3 element reports exactly the same strain
[exx, eyy, exy] = [a, d, b + c]

and stress = D(E, nu) * [exx, eyy, exy].

A constant-strain-triangle (T3/CST) element is exact for any affine
displacement field, *regardless of element size, shape, or position*. So
this check is independent of mesh geometry and directly catches any
element-size-dependent scaling error in the strain-displacement matrix
(e.g. reusing a B-matrix meant for a different, area-scaled assembly
instead of the plain kinematic strain-displacement matrix): with such a
bug, elements of different area would report different strain values for
this same imposed field, which this script would flag as a mismatch.

It bypasses the HDF5 round-trip on purpose: it runs a minimal simulation
just far enough to build the FEM dof/coords mapping, then calls
computeStrainTensorWithDisplacement / computeStressTensorWithDisplacement
directly with the synthetic displacement vector.

Usage: python3 fem_tensor_analytical.py <mesh.msh>
Exit code 0 = all elements match the analytical values; 1 = mismatch.
"""

import sys

import numpy as np

import siconos.numerics as sn

from siconos.mechanics.collision.tools import Contactor

from siconos.io.mechanics_run import MechanicsHdf5Runner, RunnerConfig
import nonos  # noqa: F401  (registers the vnative backend)

shape_filename = sys.argv[1]

backend = "vnative"
runner_config = RunnerConfig(backend)

disk_radius = 0.1

# Young's modulus / Poisson ratio must match the material passed to
# add_object() below.
E = 1e11
NU = 0.3

# Coefficients of the imposed affine displacement field
# u(x, y) = (A*x + B*y, C*x + D*y).
A, B, C, D = 1.0e-3, 2.0e-3, -1.5e-3, 0.5e-3

with MechanicsHdf5Runner(config=runner_config) as io:
    io.add_primitive_shape("DiskR", "Disk", [disk_radius])

    io.add_object(
        "disk-1",
        [Contactor("DiskR")],
        translation=[0.8, 1.5],
        orientation=[0],
        velocity=[0, 0, 0],
        mass=1,
        inertia=1,
    )

    io.add_primitive_shape("Ground-1", "Segment", (-10, 0, 10, 0))

    io.add_shape_data_from_file("Square", shape_filename)
    io.add_object("square", [Contactor("Square")], material=(2500, E, NU))

    io.add_Newton_impact_friction_nsl("contact", mu=0.5, e=0)

options = sn.SolverOptions(sn.solver_ids.SICONOS_FRICTION_2D_NSGS)
options.iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-2
options.iparam[sn.params.SICONOS_NSGS_FREEZING_CONTACT] = 10

ok = True

with MechanicsHdf5Runner(mode="r+", config=runner_config) as io:
    # Only a couple of steps are needed to trigger the FEM dof/coords
    # mapping setup inside run(); the actual simulated motion is irrelevant
    # since the displacement field checked below is synthetic.
    io.run(
        with_timer=False,
        t0=0,
        T=0.01,
        h=0.005,
        theta=0.5001,
        Newton_max_iter=1,
        set_external_forces=None,
        solver_options=options,
        numerics_verbose=False,
        output_contact_forces=False,
        output_frequency=None,
    )

    D11 = E / (1.0 - NU * NU)
    D12 = E * NU / (1.0 - NU * NU)
    D33 = E / (2.0 * (1.0 + NU))

    exx_expected = A
    eyy_expected = D
    exy_expected = B + C
    eps_expected = np.array([exx_expected, eyy_expected, exy_expected])
    sig_expected = np.array(
        [
            D11 * exx_expected + D12 * eyy_expected,
            D12 * exx_expected + D11 * eyy_expected,
            D33 * exy_expected,
        ]
    )

    if not io._fem_dof_mappings:
        print("FAIL: no FEM dynamical system found")
        ok = False

    for ds_id, mapping in io._fem_dof_mappings.items():
        fem_ds = io._nsds.dynamicalSystem(ds_id)
        fesolid = fem_ds._fesolid

        dof_indices = mapping["dof_indices"]
        coords = mapping["coords"]
        n_vertices = mapping["n_vertices"]
        num_elements = mapping["num_elements"]

        u = np.zeros(int(dof_indices.max()) + 1)
        for i in range(n_vertices):
            x, y = coords[i, 0], coords[i, 1]
            u[dof_indices[2 * i]] = A * x + B * y
            u[dof_indices[2 * i + 1]] = C * x + D * y

        epsilon_flat = np.asarray(fesolid.computeStrainTensorWithDisplacement(u))
        sigma_flat = np.asarray(fesolid.computeStressTensorWithDisplacement(u))

        if epsilon_flat.size != 3 * num_elements or sigma_flat.size != 3 * num_elements:
            print(
                f"FAIL element count mismatch for ds_id={ds_id}: "
                f"num_elements={num_elements}, "
                f"strain entries={epsilon_flat.size}, stress entries={sigma_flat.size}"
            )
            ok = False
            continue

        epsilon = epsilon_flat.reshape(num_elements, 3)
        sigma = sigma_flat.reshape(num_elements, 3)

        if not np.allclose(epsilon, eps_expected, rtol=1e-8, atol=1e-12):
            print(f"FAIL strain mismatch for ds_id={ds_id}:")
            print(f"  computed (per element):\n{epsilon}")
            print(f"  expected (constant):    {eps_expected}")
            ok = False

        if not np.allclose(sigma, sig_expected, rtol=1e-6, atol=1e-3):
            print(f"FAIL stress mismatch for ds_id={ds_id}:")
            print(f"  computed (per element):\n{sigma}")
            print(f"  expected (constant):    {sig_expected}")
            ok = False

if ok:
    print("OK: strain/stress tensors match analytical affine-field values "
          "for all elements")
    sys.exit(0)
else:
    sys.exit(1)
