"""Formulation for string-like structures,
derivation of a siconos Lagrangian DS.

Based on JSV Issanchou paper.
"""
import math
import numpy as np
import siconos.kernel as sk
import scipy.sparse as scs


class StringDS(sk.LagrangianLinearTIDS):
    """Build a string as a LagrangianLinearTIDS
    with a 'modal' spatial discretisation.
    """
    __damping_parameters_names = ['nu_air', 'rho_air', 'delta_ve', '1/qte']
    CLAMPED = 1

    def __init__(self, ndof, stiffness_coeff, density, diameter,
                 length, q0, v0, tension, damping_parameters):
        """Build string

        Parameters
        ----------
        ndof : int
            system size
        stiffness_coeff : double
            'B' parameter in (10), i.e. the 'material' part
            in eigenfrequencies. Depends on Young modulus,
            moment of inertia ...
        density : double
            linear mass density
        diameter, length, tension : double
        q0, v0 : np.ndarrays
            initial values for positions and velocities
        damping_parameters: dictionnary
            constant parameters related to damping.
            keys must be ['nu_air', 'rho_air', 'delta_ve', '1/qte']

        """
        self.ndof = ndof
        self._N = ndof - 1
        # mass = sk.SimpleMatrix(ndof, ndof, sk.SPARSE, ndof)

        # for i in xrange(ndof):
        #     mass.setValue(i, i, density)
        mass = np.identity(ndof)
        mass *= density
        self.stiffness_coeff = stiffness_coeff
        self.length = length
        self.diameter = diameter
        self.density = density
        self.c0 = math.sqrt(tension / density)
        self.tension = tension
        self.damping_parameters = damping_parameters
        msg = 'StringDS : missing parameter value for damping.'
        for name in self.__damping_parameters_names:
            assert name in self.damping_parameters.keys(), msg
        stiffness_mat, damping_mat = self.compute_linear_coeff()
        super(StringDS, self).__init__(q0, v0, mass)
        self.setCPtr(damping_mat)
        self.setKPtr(stiffness_mat)
        self.apply_boundary_conditions()

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
        indices = np.asarray(np.arange(self.ndof), dtype=np.float64)
        row = indices.reshape(1, self.ndof)
        s_mat = row * row.T
        s_mat *= (math.pi / self._N)
        s_mat = math.sqrt(2. / self.length) * np.sin(s_mat)
        # --- Omega^2 matrix ---
        omega = np.asarray([1. + self.stiffness_coeff * j ** 2
                            for j in xrange(self.ndof)], dtype=np.float64)
        coeff = (2 * math.pi * self.c0) ** 2
        omega *= coeff
        omega_mat = scs.csr_matrix((omega, (indices, indices)),
                                   shape=(self.ndof, self.ndof))
        # S.Omega^2.S-1
        stiffness_mat = s_mat * omega_mat * s_mat.T
        coeff = self.density * self.length / self._N
        stiffness_mat *= coeff

        # 2.S.Gamma.S-1
        sigma = self.compute_damping(omega / (2. * math.pi))
        sigma_mat = scs.csr_matrix((sigma, (indices, indices)),
                                   shape=(self.ndof, self.ndof))
        damping_mat = s_mat * sigma_mat * s_mat.T
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
        r = np.sqrt(math.pi * nu_air * rho_air * nu)
        r *= self.diameter
        r += nu_air
        r *= 2. * math.pi
        q_air_inv = r / (2. * math.pi * self.density * nu)
        q_ve_inv = 4 * self.length ** 2 * self.stiffness_coeff * self.density
        q_ve_inv *= delta_ve
        q_ve_inv *= nu ** 2
        q_ve_inv /= self.tension
        quality_factor_inv = q_ve_inv + q_air_inv
        quality_factor_inv += self.damping_parameters['1/qte']
        return quality_factor_inv * math.pi * nu

        # def compute_k(self):
    #     """Compute 'K' matrix = density*L/ndof * Omega**2 S^T
    #     """

    #     k_matrix = scs.csr_matrix()k.SimpleMatrix(self.ndof, self.ndof, sk.SPARSE, self.ndof)
    #     coeff = self.density * self.length / self.ndof
        
    #     it = xrange(1, self.ndof)
    #     for i, j in it, it:
    #         val = math.sqrt(2. / self.length) * math.sin(j * i * coeff)
    #         s_matrix.setValue(i, j, val)
    #     return s_matrix


    
    def computeMass(self, position=None):
        """Update inertia matrix"""
        print self._mass.nnz()
        
# class ModalString(String):
#     """String with a modal space discretisation (see Issanchou ref)
#     """

#     def __init__(self, nb_modes, **kwds):
#         """
#         """
#         # number of elements in the mesh
#         super(ModalString, self).__init__(**kwds)

#         self.nb_elements = nb_elements
#         nb_nodes = self.nb_elements + 1
#         self.stress = np.zeros((nb_elements, ), dtype=np.float64)
#         self.stress_right = np.zeros_like(self.stress)
#         self.stress_left = np.zeros_like(self.stress)
#         self.velocity = np.zeros_like(self.stress)
#         self.displacement = np.zeros((nb_nodes, ), dtype=np.float64)
#         self.velocity_right = np.zeros_like(self.stress)
#         self.velocity_left = np.zeros_like(self.stress)
#         self.celerity = math.sqrt(young_modulus / density)
#         self.rho_c = density * self.celerity
#         self.inv_rho_c = 1. / self.rho_c
#         self.inv_cross_section = 1. / cross_section
#         self.external_forces = np.zeros((nb_nodes, ), dtype=np.float64)
#         self.space_step = length / nb_elements
#         self.young_modulus = young_modulus
#         self.length = length

#     def apply_initial_state(self, u0, v0):
#         """Apply initial condition for velocity and stress
#         """
#         assert v0.size == u0.size - 1
#         assert v0.size == self.nb_elements
#         self.velocity[...] = v0
#         self.displacement[...] = u0
#         self.stress[...] = u0[1:self.nb_elements + 1] - u0[:self.nb_elements]
#         self.stress *= self.young_modulus / self.space_step

#     def compute_external_forces(self, fleft, fright, body_force):
#         """set external force distribution
#         """
#         self.external_forces[...] = body_force * self.space_step
#         self.external_forces[0] *= 0.5
#         self.external_forces[0] += fleft
#         self.external_forces[-1] *= 0.5
#         self.external_forces[-1] += fright

#     def compute_stressed_inertial_state(self, dt):
#         """Apply formula 1.17 to compute the new stressed-inertial state
#         i.e stress and velocity of each element at the end of the time step.
#         """
#         # warning : in place computation
#         self.stress *= -1
#         self.stress += self.stress_left
#         self.stress += self.stress_right
#         self.velocity *= -1
#         self.velocity += self.velocity_left
#         self.velocity += self.velocity_right
#         self.displacement += self.velocity * dt

#     def compute_inner_nodes_right_values(self):
#         """Apply formula 1.16 to compute right values of stress and velocity for inner
#         nodes
#         """
#         for j in xrange(1, self.nb_elements):
#             self.stress_right[j - 1] = self.velocity[j]
#             self.stress_right[j - 1] -= self.velocity[j - 1]
#             self.stress_right[j - 1] *= self.rho_c
#             self.stress_right[j - 1] += self.stress[j]
#             self.stress_right[j - 1] += self.stress[j - 1]
#             self.stress_right[j - 1] += self.external_forces[j] * \
#                 self.inv_cross_section
#             self.stress_right[j - 1] *= 0.5
#             self.velocity_right[j - 1] = (self.stress_right[j - 1] -
#                                           self.stress[j]) * self.rho_c +\
#                 self.velocity[j]

#             # self.velocity_right[j - 1] = self.stress[j]
#             # self.velocity_right[j - 1] -= self.stress[j - 1]
#             # self.velocity_right[j - 1] += self.external_forces[j] * \
#             #     self.inv_cross_section
#             # self.velocity_right[j - 1] *= self.inv_rho_c
#             # self.velocity_right[j - 1] += self.velocity[j]
#             # self.velocity_right[j - 1] += self.velocity[j - 1]
#             # self.velocity_right[j - 1] *= 0.5

#     def apply_continuity_and_equilibrium(self):
#         """Apply 1.15 to compute left values of stress and velocity in each element
#         """
#         self.stress_left[1:] = self.stress_right[: self.nb_elements - 1]
#         self.stress_left[1:] -= self.external_forces[1:-1] * \
#             self.inv_cross_section
#         self.velocity_left[1:] = self.velocity_right[: self.nb_elements - 1]

#     def apply_boundary_conditions(self):
#         """Apply 1.20 to compute values of stress and velocity:
#           - on the left for first element
#           - on the right for last element
#         """
#         self.stress_left[0] = self.external_forces[0] * self.inv_cross_section
#         self.stress_right[-1] = self.external_forces[-1] * \
#             self.inv_cross_section
#         self.velocity_left[0] = self.velocity[0]
#         self.velocity_left[0] -= (self.stress_left[0] -
#                                   self.stress[0]) * self.inv_rho_c
#         self.velocity_right[-1] = self.velocity[-1]
#         self.velocity_right[-1] += (self.stress_right[-1] -
#                                     self.stress[-1]) * self.inv_rho_c

#     def advance_step(self, dt):
#         """compute rod state for next time.
#         """
#         assert dt == self.space_step / self.celerity
#         self.compute_inner_nodes_right_values()
#         self.apply_continuity_and_equilibrium()
#         self.compute_stressed_inertial_state(dt)

#     def plot_state(self):
#         """plot current stress and velocity in the rod
#         """
#         x = np.linspace(0., self.length, self.nb_elements)
#         plt.subplot(2, 1, 1)
#         plt.plot(x, self.stress, label='stress')
#         plt.legend()
#         plt.subplot(2, 1, 2)
#         plt.plot(x, self.velocity, label='velocity')
#         plt.legend()
#         plt.show()

#     def __str__(self):
#         res_v = np.concatenate((self.velocity_left, self.velocity,
#                                 self.velocity_right),
#                                axis=0).reshape(3, self.nb_elements)
#         res_s = np.concatenate((self.stress_left, self.stress,
#                                 self.stress_right),
#                                axis=0).reshape(3, self.nb_elements)
#         disp = 'velocities: \n' + str(res_v) + '\n stresses:\n' + str(res_s)
#         disp += '\n displacements:\n' + str(self.displacement)
#         return disp
