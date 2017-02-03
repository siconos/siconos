"""Formulation for string-like structures,
derivation of a siconos Lagrangian DS.

Based on JSV Issanchou paper.
"""
import math
import numpy as np
import siconos.kernel as sk
# import scipy.sparse as scs
import numpywrappers as npw


class StringDS(sk.LagrangianLinearTIDS):
    """Build a string as a LagrangianLinearTIDS
    """
    __damping_parameters_names = ['nu_air', 'rho_air', 'delta_ve', '1/qte']
    CLAMPED = 1

    def __init__(self, ndof, geometry_and_material,
                 umax, imax, damping_parameters,
                 use_sparse=True, modal_form=True):
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
        use_sparse: bool, optional
            true to use sparse storage for K,C ... matrices, default=True
        modal_form: bool, optional
            true to use modal formalisation of the dynamics, default=True

        """
        self.ndof = ndof
        self._N = ndof - 1
        self.density = geometry_and_material['density']
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
        self.use_sparse = use_sparse
        self.modal_form = modal_form
        if not modal_form:
            self.use_sparse = False
        if self.modal_form:
            if use_sparse:
                self.compute_linear_coeff = self._compute_linear_coeff_sparse
            else:
                self.compute_linear_coeff = self._compute_linear_coeff_dense
            self.compute_initial_state = self._compute_initial_state_modal
        else:
            self.compute_linear_coeff = self._compute_linear_coeff_std
            self.compute_initial_state = self._compute_initial_state_std

        self.s_mat = self.compute_s_mat()
        mass, stiffness_mat, damping_mat = self.compute_linear_coeff()
        q0 = self.compute_initial_state(imax, umax)
        v0 = npw.zeros_like(q0)
        super(StringDS, self).__init__(q0, v0, mass,
                                       stiffness_mat, damping_mat)
        if not self.modal_form:
            self.apply_boundary_conditions()

    def _compute_initial_state_std(self, imax, umax):
        """Set initial positions of the string,
        assuming a triangular shape, with u[imax] = umax.
        """
        dx = self.length / self._N
        slope = umax / (imax * dx)
        q0 = [slope * i * dx for i in xrange(imax)]
        slope = umax / (imax * dx - self.length)
        q0 += [slope * (i * dx - self.length) for i in xrange(imax, self.ndof)]
        return npw.asrealarray(q0)

    def _compute_initial_state_modal(self, imax, umax):
        """Set initial positions of the string,
        assuming a triangular shape, with u[imax] = umax
        and modal form.
        """
        q0 = self._compute_initial_state_std(imax, umax)
        q0[...] = np.dot(self.s_mat.T, q0)
        coeff = self.length / self._N
        q0 *= coeff
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

    def _compute_linear_coeff_sparse(self):
        """Compute M, K and C matrices of the Lagrangian DS
        K = Omega^2
        C = 2.Gamma
        """
        mass = sk.SimpleMatrix(self.ndof, self.ndof, sk.SPARSE, self.ndof)
        for i in xrange(self.ndof):
            mass.setValue(i, i, 1.)
        indices = np.arange(self.ndof)
        # --- Omega^2 matrix ---
        omega = npw.asrealarray([1. + self.stiffness_coeff * j ** 2
                                 for j in xrange(self.ndof)])
        coeff = (indices * math.pi * self.c0 / self.length) ** 2
        omega *= coeff
        stiffness_mat = sk.SimpleMatrix(self.ndof, self.ndof,
                                        sk.SPARSE, self.ndof)
        for i in xrange(self.ndof):
            stiffness_mat.setValue(i, i, omega[i])

        # 2.S.Gamma.S-1
        sigma = self.compute_damping(np.sqrt(omega) / (2. * math.pi))
        damping_mat = sk.SimpleMatrix(self.ndof, self.ndof,
                                      sk.SPARSE, self.ndof)
        for i in xrange(self.ndof):
            damping_mat.setValue(i, i, 2. * sigma[i])
        return mass, stiffness_mat, damping_mat

    def _compute_linear_coeff_dense(self):
        """Compute M, K and C matrices of the Lagrangian DS
        K = Omega^2
        C = 2.Gamma
        """
        mass = np.identity(self.ndof, dtype=np.float64)

        # --- Omega^2 matrix ---
        omega = npw.asrealarray([1. + self.stiffness_coeff * j ** 2
                                 for j in xrange(self.ndof)])
        indices = np.arange(self.ndof)
        coeff = (indices * math.pi * self.c0 / self.length) ** 2
        omega *= coeff
        stiffness_mat = npw.asrealarray(np.diag(omega))

        # 2.S.Gamma.S-1
        sigma = self.compute_damping(np.sqrt(omega) / (2. * math.pi))
        damping_mat = npw.asrealarray(np.diag(sigma))
        damping_mat *= 2.
        return mass, stiffness_mat, damping_mat

    def _compute_linear_coeff_std(self):
        """Compute K and C matrices of the Lagrangian DS
        K = S.Omega^2.S^-1
        C = 2S.Gamma.S^-1
        """
        mass, omega_mat, sigma_mat = self._compute_linear_coeff_dense()
        # S.Omega^2.S-1
        stiffness_mat = np.dot(self.s_mat, np.dot(omega_mat, self.s_mat.T))
        coeff = self.length / self._N
        stiffness_mat *= coeff

        # 2.S.Gamma.S-1
        damping_mat = np.dot(self.s_mat, np.dot(sigma_mat, self.s_mat.T))
        damping_mat *= coeff
        return mass, stiffness_mat, damping_mat

    def compute_s_mat(self):
        """Compute 'S' matrix (normal modes)
        """
        indices = np.arange(self.ndof)
        row = indices.reshape(1, self.ndof)
        s_mat = npw.asrealarray(row * row.T)
        s_mat *= (math.pi / self._N)
        s_mat[...] = math.sqrt(2. / self.length) * np.sin(s_mat)
        return s_mat

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
