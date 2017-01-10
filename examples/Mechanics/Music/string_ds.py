"""Formulation for string-like structures,
derivation of a siconos Lagrangian DS.

Based on JSV Issanchou paper.
"""
import math
import numpy as np
import siconos.kernel as sk
import scipy.sparse as scs
import numpywrappers as npw


class StringDS(sk.LagrangianLinearTIDS):
    """Build a string as a LagrangianLinearTIDS
    with a 'modal' spatial discretisation.
    """
    __damping_parameters_names = ['nu_air', 'rho_air', 'delta_ve', '1/qte']
    CLAMPED = 1

    def __init__(self, ndof, geometry_and_material,
                 umax, imax, damping_parameters):
        """Build string

        Parameters
        ----------
        ndof : int
            system size
        geometry_and_material : dictionnary
            description of the string, keys must be:
            ['length', 'diameter', 'B', 'tension', 'density']
        umax : double
            initial value of string maximum position
            (assuming triangular shape)
        imax : int
            index of maximum initial position
            (i.e. u[imax, t=0] = umax)
        damping_parameters: dictionnary
            constant parameters related to damping.
            keys must be ['nu_air', 'rho_air', 'delta_ve', '1/qte']

        """
        self.ndof = ndof
        self._N = ndof - 1
        # mass = sk.SimpleMatrix(ndof, ndof, sk.SPARSE, ndof)

        # for i in xrange(ndof):
        #     mass.setValue(i, i, density)
        mass = np.identity(ndof, dtype=np.float64)
        self.density = geometry_and_material['density']
        #mass *= self.density
        self.stiffness_coeff = geometry_and_material['B']
        # B is a function of Young Modulus (E),
        # moment of inertia ...
        # B = pi^2EI / (TL^2)
        self.length = geometry_and_material['length']
        self.diameter = geometry_and_material['diameter']
        self.tension = geometry_and_material['tension']
        self.c0 = math.sqrt(self.tension / self.density)
        self.damping_parameters = damping_parameters
        msg = 'StringDS : missing parameter value for damping.'
        for name in self.__damping_parameters_names:
            assert name in self.damping_parameters.keys(), msg
        stiffness_mat, damping_mat = self.compute_linear_coeff()
        q0 = self.compute_initial_state(imax, umax)
        v0 = npw.zeros_like(q0)
        super(StringDS, self).__init__(q0, v0, mass)
        self.setCPtr(damping_mat)
        self.setKPtr(stiffness_mat)
        self.apply_boundary_conditions()

    def compute_initial_state(self, imax, umax):
        """Set initial positions of the string,
        assuming a triangular shape, with u[imax] = umax.
        """
        dx = self.length / self._N
        slope = umax / (imax * dx)
        q0 = [slope * i * dx for i in xrange(imax)]
        slope = umax / (imax * dx - self.length)
        q0 += [slope * (i * dx - self.length) for i in xrange(imax, self.ndof)]
        return npw.asrealarray(q0)

    def eigenfreq(self, j):
        """Compute eigenfrequency number j
        """
        return 0.5 * j * self.c0 / self.length * \
            math.sqrt(1 + self.stiffness_coeff * j)

    def apply_boundary_conditions(self, bc_type=None):
        """Create and apply boundary conditions
        """
        if bc_type is self.CLAMPED or bc_type is None:
            bc_clamped_indices = [0, self.ndof - 1]
            # For each index in bc_clamped_indices
            # enforce velocity_index = 0.
            boundaries = sk.FixedBC(bc_clamped_indices)
            self.setBoundaryConditions(boundaries)
        else:
            raise AttributeError('Unknown boundary type')

    def compute_linear_coeff(self):
        """Compute K and C matrices of the Lagrangian DS
        K = S.Omega^2.S^-1
        C = 2S.Gamma.S^-1
        """
        # --- Compute 'S' matrix (normal modes) ---
        indices = np.arange(self.ndof)
        row = indices.reshape(1, self.ndof)
        s_mat = npw.asrealarray(row * row.T)
        s_mat *= (math.pi / self._N)
        s_mat = math.sqrt(2. / self.length) * np.sin(s_mat)
        # --- Omega^2 matrix ---
        omega = npw.asrealarray([1. + self.stiffness_coeff * j ** 2
                                 for j in xrange(self.ndof)])
        coeff = (indices * math.pi * self.c0 / self.length) ** 2
        omega *= coeff
        # omega_mat = scs.csr_matrix((omega, (indices, indices)),
        #                            shape=(self.ndof, self.ndof))
        omega_mat = np.diag(omega)
        #omega[0] = 1.
        # S.Omega^2.S-1
        stiffness_mat = np.dot(s_mat, np.dot(omega_mat, s_mat.T))
        coeff = self.length / self._N
        stiffness_mat *= coeff

        # 2.S.Gamma.S-1
        sigma = self.compute_damping(np.sqrt(omega) / (2. * math.pi))
        # sigma_mat = scs.csr_matrix((sigma, (indices, indices)),
        #                            shape=(self.ndof, self.ndof))
        sigma_mat = np.diag(sigma)
        damping_mat = np.dot(s_mat, np.dot(sigma_mat, s_mat.T))
        coeff = 2. * coeff
        damping_mat *= coeff
        return stiffness_mat, damping_mat

    def compute_damping(self, nu):
        """Compute inv of quality factor

        Parameters
        ----------
        nu : 1D numpy array
            vector of eigenfrequencies

        Returns
        -------
        vector of sigma_j values (i.e. diagonal of gamma in (11))
        """
        # Compute quality factor as defined in
        # eq (12)
        nu_air = self.damping_parameters['nu_air']
        rho_air = self.damping_parameters['rho_air']
        delta_ve = self.damping_parameters['delta_ve']
        r = np.sqrt(math.pi * nu_air * rho_air * nu[1:])
        r *= self.diameter
        r += nu_air
        r *= 2. * math.pi
        q_air_inv = r / (2. * math.pi * self.density * nu[1:])
        q_ve_inv = 4 * self.length ** 2 * self.stiffness_coeff * self.density
        q_ve_inv *= delta_ve
        q_ve_inv *= nu[1:] ** 2
        q_ve_inv /= self.tension
        quality_factor_inv = q_ve_inv + q_air_inv
        quality_factor_inv += self.damping_parameters['1/qte']
        return np.insert(quality_factor_inv * math.pi * nu[1:], 0, 0)
