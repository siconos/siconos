# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
"""Plasticity problems (numerics/src/Plasticity, docs/devel/PLASTICITY.md).

Solves with plasticity_2d_driver:
- Drucker-Prager: analytic solutions with M = I (cohesion, non-trivial case with
  q_N >= 0), coupled problems with dilatancy, all the local solvers;
- Von Mises plane strain (cylinder): pure shear and plane strain uniaxial loading;
- Von Mises plane stress (ellipsoid, solved as a ball): manufactured solutions with a
  non-symmetric M, Voigt / Mandel / cone coordinates, all the local solvers;
- Cartesian coordinates: same solution as the cone problem.
The maps between coordinates and the plane stress formula are written here from their
definitions, independently of the library.
"""

import numpy as np
import pytest
import siconos.numerics as sn

sn.numerics_set_verbose(0)

TOL = sn.params.SICONOS_DPARAM_enum.SICONOS_DPARAM_TOL
RESIDU = sn.params.SICONOS_DPARAM_enum.SICONOS_DPARAM_RESIDU
MAX_ITER = sn.params.SICONOS_IPARAM_enum.SICONOS_IPARAM_MAX_ITER
IDS = sn.solver_ids.PLASTICITY_SOLVER
COORDS = sn.PlasticityStressCoordinates
HYP = sn.PlasticityVonMisesHypothesis

DP_LOCAL_SOLVERS = [
    None,  # default of the NSGS (semi-smooth Newton)
    IDS.PLASTICITY_2D_ONECONE_NSN,
    IDS.PLASTICITY_2D_ONECONE_ProjectionOnCone,
    IDS.PLASTICITY_2D_ONECONE_ProjectionOnConeWithLocalIteration,
    IDS.PLASTICITY_2D_ONECONE_DRUCKERPRAGER_RETURN_MAPPING,
    IDS.PLASTICITY_2D_ONECONE_DRUCKERPRAGER_NSN,
]
VM_LOCAL_SOLVERS = [
    None,
    IDS.PLASTICITY_2D_ONECONE_VONMISES_RETURN_MAPPING,
    IDS.PLASTICITY_2D_ONECONE_VONMISES_FP,
    IDS.PLASTICITY_2D_ONECONE_VONMISES_NSN,
    IDS.PLASTICITY_2D_ONECONE_ProjectionOnCone,
]


def solver_name(local):
    return "default" if local is None else local.name


def solve(problem, local=None, tol=1e-14, itermax=20000):
    n = 3 * problem.numberOfCones
    stress, strainrate = np.zeros(n), np.zeros(n)
    options = sn.SolverOptions(int(IDS.PLASTICITY_2D_NSGS))
    options.dparam[TOL] = tol
    options.iparam[MAX_ITER] = itermax
    if local is not None:
        options.update_internal(0, int(local))
    info = sn.plasticity_2d_driver(problem, stress, strainrate, options)
    return info, stress, strainrate, options.dparam[RESIDU]


def error(problem, stress):
    w = np.zeros_like(stress)
    _, err = sn.plasticity_2d_compute_error(problem, stress.copy(), w, 1e-12)
    return err


def random_coupled_matrix(nc, rng):
    """Non-symmetric matrix with a positive definite symmetric part, coupling neighbours."""
    n = 3 * nc
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if abs(i // 3 - j // 3) <= 1:
                M[i, j] = (0.3 if i // 3 == j // 3 else 0.1) * rng.uniform(-1, 1)
    return M + 2.0 * np.eye(n)


# --------------------------------------------------------------------------------------
# Coordinates (definitions)
# --------------------------------------------------------------------------------------
H = np.sqrt(0.5)
# sigma_c = P s for one element, Mandel s = (sxx, syy, sqrt2 sxy), Voigt s = (sxx, syy, sxy)
P_MANDEL = np.array([[-H, -H, 0.0], [H, -H, 0.0], [0.0, 0.0, 1.0]])
P_VOIGT = P_MANDEL @ np.diag([1.0, 1.0, np.sqrt(2.0)])


def block_diag(A, nc):
    return np.kron(np.eye(nc), A)


def stress_vector(coords, c):
    """Stress vector of (sxx, syy, sxy) in the given coordinates."""
    sxx, syy, sxy = c
    if coords == COORDS.CARTESIAN_VOIGT:
        return np.array([sxx, syy, sxy])
    if coords == COORDS.CARTESIAN_MANDEL:
        return np.array([sxx, syy, np.sqrt(2.0) * sxy])
    return P_MANDEL @ np.array([sxx, syy, np.sqrt(2.0) * sxy])


def cartesian(coords, s):
    """(sxx, syy, sxy) of a stress vector."""
    if coords == COORDS.CARTESIAN_VOIGT:
        return s.copy()
    if coords == COORDS.CARTESIAN_MANDEL:
        return np.array([s[0], s[1], s[2] / np.sqrt(2.0)])
    m = P_MANDEL.T @ s
    return np.array([m[0], m[1], m[2] / np.sqrt(2.0)])


def plane_stress_gradient(coords, s):
    """Gradient of f^2/2, f^2 = sxx^2 + syy^2 - sxx syy + 3 sxy^2, w.r.t. the vector s."""
    if coords == COORDS.CARTESIAN_VOIGT:
        return np.array([s[0] - 0.5 * s[1], s[1] - 0.5 * s[0], 3.0 * s[2]])
    if coords == COORDS.CARTESIAN_MANDEL:
        return np.array([s[0] - 0.5 * s[1], s[1] - 0.5 * s[0], 1.5 * s[2]])
    return np.array([0.5 * s[0], 1.5 * s[1], 1.5 * s[2]])  # sN^2/2 + 3/2 |sT|^2


def plane_stress_equivalent(c):
    sxx, syy, sxy = c
    return np.sqrt(sxx**2 + syy**2 - sxx * syy + 3.0 * sxy**2)


# --------------------------------------------------------------------------------------
# Problem creation
# --------------------------------------------------------------------------------------
M1 = np.eye(3)
q1 = np.array([-1.0, 1.0, 3.0])
eta1 = np.array([0.1])
theta1 = np.array([0.1])


def test_plasticity_problem_creation():
    P2D = sn.PlasticityProblem_new(3, M1, q1, eta1, theta1)
    assert P2D.dimension == 3
    assert P2D.numberOfCones == 1
    assert P2D.model_type == sn.PlasticityModelType.DRUCKER_PRAGER
    assert P2D.stress_coordinates == COORDS.CONE


def test_plasticity_problem_2d_alias():
    P2D = sn.Plasticity2DProblem_new(3, M1, q1, eta1, theta1)
    assert P2D.dimension == 3
    assert P2D.numberOfCones == 1


def test_solver_ids():
    assert sn.PLASTICITY_2D_NSGS == 20000
    assert sn.MOHR_COULOMB_2D_NSGS == 20000  # backward compatibility
    assert sn.SICONOS_MC2D_NSGS == 20000  # backward compatibility
    assert sn.PLASTICITY_2D_ONECONE_NSN == 20050


def test_solver_options_create():
    options = sn.SolverOptions(sn.solver_ids.PLASTICITY_2D_NSGS)
    assert options is not None


def test_invalid_sizes():
    with pytest.raises(ValueError):
        sn.PlasticityProblem_new(3, np.eye(4), q1, eta1, theta1)
    with pytest.raises(ValueError):
        sn.PlasticityProblem_new(3, M1, q1, eta1, np.array([0.1, 0.2]))
    problem = sn.PlasticityProblem_new(3, M1, q1, eta1, theta1)
    with pytest.raises(ValueError):
        sn.plasticity_2d_driver(problem, np.zeros(2), np.zeros(3))


def test_matrix_layout():
    """M is used as given (row i of the numpy array = row i of M), C or Fortran order."""
    rng = np.random.default_rng(0)
    M = random_coupled_matrix(2, rng)
    q = rng.uniform(-1, 1, 6)
    s = rng.uniform(-1, 1, 6)
    for Mi in (M, np.asfortranarray(M)):
        problem = sn.PlasticityProblem_new_VonMises(3, Mi, q, np.array([1.0, 1.0]))
        w = np.zeros(6)
        sn.plasticity_2d_compute_error(problem, s.copy(), w, 1e-12)
        assert np.allclose(w, M @ s + q, rtol=0, atol=1e-14)


def test_error_scaling():
    """Scaling of the error as in fc3d: the MIXED scaling makes it relative, homogeneous of
    degree 0 in (q, stress)."""
    rng = np.random.default_rng(4)
    M = random_coupled_matrix(3, rng)
    q = rng.uniform(-1, 1, 9)
    stress = rng.uniform(4, 6, 9)
    P = sn.params
    mixed = sn.SolverOptions(int(IDS.PLASTICITY_2D_NSGS))
    mixed.iparam[P.PLASTICITY_SCALING_REDUCTION] = int(
        P.PLASTICITY_SCALING_REDUCTION_MIXED
    )
    absolute = sn.SolverOptions(int(IDS.PLASTICITY_2D_NSGS))
    absolute.iparam[P.PLASTICITY_SCALING_NORM] = int(P.PLASTICITY_SCALING_NO)
    errors = {}
    for lam in (1.0, 1e4):
        problem = sn.PlasticityProblem_new(3, M, lam * q, np.full(3, 0.4), np.zeros(3))
        for name, opt in (("mixed", mixed), ("absolute", absolute)):
            w = np.zeros(9)
            _, errors[name, lam] = sn.plasticity_2d_compute_error(
                problem, lam * stress, w, 1e-12, opt
            )
    assert abs(errors["absolute", 1e4] / errors["absolute", 1.0] - 1e4) < 1e-6
    assert abs(errors["mixed", 1e4] / errors["mixed", 1.0] - 1.0) < 1e-10
    # rho_e = 1/||M_ee||: the error does not depend on the units of the strain rate
    rho = sn.SolverOptions(int(IDS.PLASTICITY_2D_NSGS))
    rho.iparam[P.PLASTICITY_SCALING_RHO] = int(P.PLASTICITY_SCALING_RHO_ELEMENT)
    for kappa in (1.0, 1e-6):
        problem = sn.PlasticityProblem_new(
            3, kappa * M, kappa * q, np.full(3, 0.4), np.zeros(3)
        )
        _, errors["rho", kappa] = sn.plasticity_2d_compute_error(
            problem, stress.copy(), np.zeros(9), 1e-12, rho
        )
    assert abs(errors["rho", 1e-6] / errors["rho", 1.0] - 1.0) < 1e-10


def test_default_options():
    problem = sn.PlasticityProblem_new(3, M1, q1, eta1, theta1)
    stress, w = np.zeros(3), np.zeros(3)
    info = sn.plasticity_2d_driver(problem, stress, w)
    assert info == 0


# --------------------------------------------------------------------------------------
# Drucker-Prager
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("local", DP_LOCAL_SOLVERS, ids=solver_name)
def test_drucker_prager_cohesion_identity(local):
    """M = I, eta = 0.5, theta = 0, c = 0.5: sigma = P_K(-q) (tension within the
    cohesion, apex, sliding), see plasticity_2d_cohesion_test.c."""
    q = np.array([0.1, 0.0, 0.0, 2.0, 0.0, 0.0, -1.0, -2.0, 0.0])
    expected = np.array([-0.1, 0.0, 0.0, -1.0, 0.0, 0.0, 1.4, 1.2, 0.0])
    eta, theta, c = np.full(3, 0.5), np.zeros(3), np.full(3, 0.5)
    problem = sn.PlasticityProblem_new(3, np.eye(9), q, eta, theta, cohesion=c)
    info, stress, w, _ = solve(problem, local)
    assert info == 0
    assert np.allclose(stress, expected, rtol=0, atol=1e-13)
    assert np.allclose(w, stress + q, rtol=0, atol=1e-13)
    assert error(problem, stress) < 1e-12


@pytest.mark.parametrize("local", DP_LOCAL_SOLVERS, ids=solver_name)
def test_drucker_prager_not_trivial(local):
    """q_N >= 0 but sigma = 0 is not a solution (theta < eta): sigma = P_K(-q)."""
    problem = sn.PlasticityProblem_new(
        3, np.eye(3), np.array([0.1, 1.0, 0.0]), np.array([0.5]), np.array([0.0])
    )
    info, stress, _, _ = solve(problem, local)
    assert info == 0
    assert np.allclose(stress, [0.32, -0.16, 0.0], rtol=0, atol=1e-13)


def test_drucker_prager_dilatancy_coupled():
    """Coupled non-symmetric problem with dilatancy and cohesion: all the local solvers
    converge to the same solution, with a small error."""
    rng = np.random.default_rng(3)
    nc = 5
    M = random_coupled_matrix(nc, rng)
    q = rng.uniform(-1, 1, 3 * nc)
    eta = 0.3 + 0.05 * np.arange(nc)
    theta = 0.1 * (np.arange(nc) % 2)
    c = 0.05 * np.arange(nc)
    solutions = []
    for local in DP_LOCAL_SOLVERS:
        problem = sn.PlasticityProblem_new(3, M, q, eta, theta, cohesion=c)
        info, stress, _, _ = solve(problem, local)
        assert info == 0, solver_name(local)
        assert error(problem, stress) < 1e-12
        solutions.append(stress)
    # some cones on the yield surface |sigma_T| = eta sigma_N + c
    s = solutions[0].reshape(nc, 3)
    on_surface = np.abs(np.hypot(s[:, 1], s[:, 2]) - (eta * s[:, 0] + c)) < 1e-10
    assert on_surface.any()
    for stress in solutions[1:]:
        assert np.allclose(stress, solutions[0], rtol=0, atol=1e-10)


# --------------------------------------------------------------------------------------
# Von Mises, plane strain (cylinder)
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("local", VM_LOCAL_SOLVERS, ids=solver_name)
@pytest.mark.parametrize("coords", [COORDS.CARTESIAN_VOIGT, COORDS.CARTESIAN_MANDEL])
def test_von_mises_plane_strain(coords, local):
    """M = I: pure shear yields at sxy = sigma_y/sqrt3, plane strain uniaxial loading at
    sxx - syy = 2 sigma_y/sqrt3."""
    sigma_y, tau = 3.0, 10.0
    q = np.array([0.0, 0.0, -tau, -tau, 0.0, 0.0])
    problem = sn.PlasticityProblem_new_VonMises(3, np.eye(6), q, np.full(2, sigma_y))
    problem.stress_coordinates = coords
    info, stress, _, _ = solve(problem, local)
    assert info == 0
    shear = cartesian(coords, stress[:3])
    uniaxial = cartesian(coords, stress[3:])
    assert np.allclose(shear, [0.0, 0.0, sigma_y / np.sqrt(3.0)], rtol=0, atol=1e-12)
    assert abs(uniaxial[0] - uniaxial[1] - 2.0 * sigma_y / np.sqrt(3.0)) < 1e-12
    assert error(problem, stress) < 1e-12


# --------------------------------------------------------------------------------------
# Von Mises, plane stress (ellipsoid, solved as a ball)
# --------------------------------------------------------------------------------------
TARGETS = [
    ([1.0, 0.0, 0.0], 0.5),  # uniaxial tension
    ([0.0, -1.0, 0.0], 0.6),  # uniaxial compression
    ([1.0, 1.0, 0.0], 0.7),  # equibiaxial tension
    ([-1.0, -1.0, 0.0], 0.8),  # equibiaxial compression
    ([0.0, 0.0, 1.0 / np.sqrt(3.0)], 0.9),  # pure shear
    ([0.3, -0.2, 0.1], 0.0),  # interior point (elastic)
]


def manufactured_plane_stress(coords, rng):
    """q = -M s - lambda grad(f^2/2)(s): the unique solution is s (M non-symmetric with a
    positive definite symmetric part)."""
    nc = len(TARGETS) + 1
    sigma_y = 1.0 + 0.5 * np.arange(nc)
    targets = [np.array(t) for t, _ in TARGETS]
    lambdas = [lam for _, lam in TARGETS]
    a, b, c = rng.uniform(-1, 1, 3)  # a random point of the yield surface
    targets.append(np.array([a, b, c]) / plane_stress_equivalent([a, b, c]))
    lambdas.append(0.4)
    s = np.concatenate(
        [stress_vector(coords, sy * t) for sy, t in zip(sigma_y, targets, strict=True)]
    )
    g = np.concatenate(
        [
            lam * plane_stress_gradient(coords, s[3 * e : 3 * e + 3])
            for e, lam in enumerate(lambdas)
        ]
    )
    M = random_coupled_matrix(nc, rng)
    q = -M @ s - g
    problem = sn.PlasticityProblem_new_VonMises(3, M, q, sigma_y, HYP.PLANE_STRESS)
    problem.stress_coordinates = coords
    return problem, s, sigma_y, lambdas


@pytest.mark.parametrize("local", VM_LOCAL_SOLVERS, ids=solver_name)
@pytest.mark.parametrize(
    "coords", [COORDS.CARTESIAN_VOIGT, COORDS.CARTESIAN_MANDEL, COORDS.CONE]
)
def test_von_mises_plane_stress(coords, local):
    problem, s_exact, sigma_y, lambdas = manufactured_plane_stress(
        coords, np.random.default_rng(11)
    )
    info, stress, _, _ = solve(problem, local)
    assert info == 0
    assert np.linalg.norm(stress - s_exact) < 1e-11 * np.linalg.norm(s_exact)
    for e, lam in enumerate(lambdas):
        if lam > 0.0:  # on the plane stress yield surface
            c = cartesian(coords, stress[3 * e : 3 * e + 3])
            assert abs(plane_stress_equivalent(c) / sigma_y[e] - 1.0) < 1e-11
    assert error(problem, stress) < 1e-11


def test_von_mises_plane_stress_differs_from_plane_strain():
    """Equibiaxial tension: on the plane stress yield surface, inside the plane strain
    cylinder (pure hydrostatic in-plane stress)."""
    q = -stress_vector(COORDS.CARTESIAN_VOIGT, [2.0, 2.0, 0.0])
    sigma_y = np.array([1.0])
    plane_stress = sn.PlasticityProblem_new_VonMises(
        3, np.eye(3), q, sigma_y, HYP.PLANE_STRESS
    )
    plane_strain = sn.PlasticityProblem_new_VonMises(3, np.eye(3), q, sigma_y)
    for problem in (plane_stress, plane_strain):
        problem.stress_coordinates = COORDS.CARTESIAN_VOIGT
    _, s_stress, _, _ = solve(plane_stress)
    _, s_strain, _, _ = solve(plane_strain)
    assert np.allclose(s_stress, [1.0, 1.0, 0.0], rtol=0, atol=1e-12)
    assert np.allclose(s_strain, [2.0, 2.0, 0.0], rtol=0, atol=1e-12)


# --------------------------------------------------------------------------------------
# Cartesian coordinates
# --------------------------------------------------------------------------------------
@pytest.mark.parametrize("coords", [COORDS.CARTESIAN_VOIGT, COORDS.CARTESIAN_MANDEL])
def test_cartesian_drucker_prager(coords):
    """The Cartesian problem M_s = P^T M_c P, q_s = P^T q_c gives s = P^{-1} sigma_c."""
    rng = np.random.default_rng(5)
    nc = 4
    Mc = random_coupled_matrix(nc, rng)
    qc = rng.uniform(-1, 1, 3 * nc)
    eta = np.full(nc, 0.4)
    theta = np.full(nc, 0.1)
    c = np.full(nc, 0.05)
    P = block_diag(P_VOIGT if coords == COORDS.CARTESIAN_VOIGT else P_MANDEL, nc)
    cone = sn.PlasticityProblem_new(3, Mc, qc, eta, theta, cohesion=c)
    cart = sn.PlasticityProblem_new(3, P.T @ Mc @ P, P.T @ qc, eta, theta, cohesion=c)
    cart.stress_coordinates = coords
    info_c, sigma_c, w_c, _ = solve(cone)
    info_s, s, w_s, _ = solve(cart)
    assert info_c == 0 and info_s == 0
    assert np.allclose(s, np.linalg.solve(P, sigma_c), rtol=0, atol=1e-12)
    assert np.allclose(w_s, P.T @ w_c, rtol=0, atol=1e-12)
    assert abs(s @ w_s - sigma_c @ w_c) < 1e-12  # power preserved
    assert error(cart, s) < 1e-12


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# --------------------------------------------------------------------------------------
# Global semi-smooth Newton (PLASTICITY_2D_NSN)
# --------------------------------------------------------------------------------------
def solve_nsn(problem, tol=1e-12, itermax=100):
    n = 3 * problem.numberOfCones
    stress, strainrate = np.zeros(n), np.zeros(n)
    options = sn.SolverOptions(int(IDS.PLASTICITY_2D_NSN))
    options.dparam[TOL] = tol
    options.iparam[MAX_ITER] = itermax
    info = sn.plasticity_2d_driver(problem, stress, strainrate, options)
    return info, stress


def test_nsn_drucker_prager_dilatancy_coupled():
    """NSN (natural-map merit: non-symmetric, dilatancy and cohesion) = NSGS."""
    rng = np.random.default_rng(3)
    nc = 5
    M = random_coupled_matrix(nc, rng)
    q = rng.uniform(-1, 1, 3 * nc)
    eta = 0.3 + 0.05 * np.arange(nc)
    theta = 0.1 * (np.arange(nc) % 2)
    c = 0.05 * np.arange(nc)
    problem = sn.PlasticityProblem_new(3, M, q, eta, theta, cohesion=c)
    info, stress = solve_nsn(problem)
    assert info == 0
    assert error(problem, stress) < 1e-12
    _, stress_nsgs, _, _ = solve(problem)
    assert np.allclose(stress, stress_nsgs, rtol=0, atol=1e-10)


@pytest.mark.parametrize("coords", [COORDS.CARTESIAN_VOIGT, COORDS.CONE])
def test_nsn_von_mises_plane_stress(coords):
    """NSN (QP merit) on the manufactured plane stress problem: exact solution."""
    problem, s_exact, _, _ = manufactured_plane_stress(
        coords, np.random.default_rng(11)
    )
    info, stress = solve_nsn(problem)
    assert info == 0
    assert np.linalg.norm(stress - s_exact) < 1e-11 * np.linalg.norm(s_exact)
    assert error(problem, stress) < 1e-11
