"""Run the same 'model' but with different ds types
and different osi types.
Ensure that:
LagrangianTIDS + MoreauJeanOSI,
LagrangianLinearDiagonalDS + MoreauJeanOSI,
LagrangianLinearDiagonalDS + MoreauJeanBilbaoOSI
"""

import siconos.kernel as sk
import numpy as np


def test_lagrangian_and_osis():
    """Tests osi/lagrangian combinations
    """
    # --- The dynamical systems ---
    # dimension
    ndof = 3
    # initial conditons
    q0 = np.zeros(ndof, dtype=np.float64)
    v0 = np.zeros_like(q0)
    q0[0] = 1.
    # mass, stiffness and damping matrices
    # (diagonal values only)
    mass_diag = np.ones(ndof, dtype=np.float64)
    k_diag = np.zeros_like(mass_diag)
    c_diag = np.zeros_like(mass_diag)
    sigma = 0.2
    omega2 = 1.2
    k_diag[...] = omega2
    c_diag[...] = sigma

    ds_list = {}
    # - LagrangianLinearTIDS -
    ds_list['LTIDS+MJ'] = sk.LagrangianLinearTIDS(q0, v0, np.diag(mass_diag),
                                                  np.diag(k_diag),
                                                  np.diag(c_diag))
    # - LagrangianLinearDiagonalDS -
    ds_list['LLDDS+MJ'] = sk.LagrangianLinearDiagonalDS(q0, v0, k_diag,
                                                        c_diag,
                                                        mass_diag)
    ds_list['LLDDS+MJB'] = sk.LagrangianLinearDiagonalDS(q0, v0, k_diag,
                                                         c_diag,
                                                         mass_diag)
    # no mass (i.e. implicitely equal to Id)
    ds_list['LLDDS+MJB2'] = sk.LagrangianLinearDiagonalDS(q0, v0,
                                                          k_diag, c_diag)

    nb_ds = len(ds_list)
    for ds in ds_list.values():
        ds.computeForces(0., q0, v0)

    assert np.allclose(ds_list['LTIDS+MJ'].forces(),
                       ds_list['LTIDS+MJ'].forces())

    # --- Interactions ---
    cor = 0.9

    nslaw = sk.NewtonImpactNSL(cor)
    hmat = np.zeros((1, ndof), dtype=np.float64)
    hmat[0, 0] = 1.
    relation = sk.LagrangianLinearTIR(hmat)

    interactions = []
    for k in range(nb_ds):
        interactions.append(sk.Interaction(nslaw, relation))

    # --- Model ---
    tinit = 0.
    tend = 3.
    model = sk.Model(tinit, tend)

    # - insert ds into the model and link them with their interaction -
    nsds = model.nonSmoothDynamicalSystem()
    ninter = 0
    for ds in ds_list.values():
        nsds.insertDynamicalSystem(ds)
        nsds.link(interactions[ninter], ds)
        ninter += 1

    # --- Simulation ---

    # - (1) OneStepIntegrators --
    theta = 0.5000001
    standard = sk.MoreauJeanOSI(theta)
    bilbao = sk.MoreauJeanBilbaoOSI()
    #-- (2) Time discretisation --
    time_step = 1e-4
    td = sk.TimeDiscretisation(tinit, time_step)
    # -- (3) one step non smooth problem
    lcp = sk.LCP()
    # -- (4) Simulation setup with (1) (2) (3)
    simu = sk.TimeStepping(td, standard, lcp)
    # extra osi must be explicitely inserted into simu and linked to ds
    simu.setNonSmoothDynamicalSystemPtr(model.nonSmoothDynamicalSystem())
    simu.prepareIntegratorForDS(bilbao, ds_list['LLDDS+MJB'], model, tinit)
    simu.prepareIntegratorForDS(bilbao, ds_list['LLDDS+MJB2'], model, tinit)

    # link simu and model, initialize
    model.setSimulation(simu)
    model.initialize()

    positions = []
    velocities = []
    for ds in ds_list.values():
        positions.append(ds.q())
        velocities.append(ds.velocity())

    # number of time steps and buffer used to save results
    nb_steps = int((tend - tinit) / time_step) + 1
    data = np.zeros((nb_steps, 1 + 2 * len(ds_list)), dtype=np.float64)
    data[0, 0] = tinit

    for k in range(nb_ds):
        data[0, k * 2 + 1] = positions[k][0]
        data[0, k * 2 + 2] = velocities[k][0]

    current_iter = 1
    # --- Time loop ---
    while simu.hasNextEvent():
        simu.computeOneStep()
        data[current_iter, 0] = simu.nextTime()
        for k in range(nb_ds):
            data[current_iter, k * 2 + 1] = positions[k][0]
            data[current_iter, k * 2 + 2] = velocities[k][0]
        simu.nextStep()
        current_iter += 1

    # -- check results
    ref_pos = data[:, 1]
    ref_vel = data[:, 2]
    for k in range(1, nb_ds):
        assert np.allclose(ref_pos, data[:, k * 2 + 1], atol=time_step)
        assert np.allclose(ref_vel, data[:, k * 2 + 2], atol=time_step)
