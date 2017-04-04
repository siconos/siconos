"""Definition of string dynamics (class StringDS)
and fret interaction (class Fret) as objects derived
from Siconos classes.
"""
import math
import numpy as np
import siconos.kernel as sk
# import scipy.sparse as scs
import numpywrappers as npw


class StringDS(sk.LagrangianLinearDiagonalDS):
    """Formulation for string-like structures,
    Based on LagrangianLinearDiagonalDS from Siconos

    Operators formulas from JSV Issanchou paper.
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


        Notes
        -----
        * DS.stiffness = diag(omega**2)
        * DS.damping  = diag(2 * sigma)
        * DS.mass = identity

        """
        self.n_modes = ndof
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
        self.s_mat = self.compute_s_mat()
        mass, stiffness_mat, damping_mat = self.compute_linear_coeff()
        q0 = self.compute_initial_state_modal(imax, umax)
        v0 = npw.zeros_like(q0)
        super(StringDS, self).__init__(q0, v0, mass,
                                       stiffness_mat, damping_mat)

    def _compute_initial_state_std(self, imax, umax):
        """Set initial positions of the string,
        assuming a triangular shape, with u[imax] = umax.
        """
        dx = self.length / (self.n_modes + 1)
        slope = umax / (imax * dx)
        q0 = [slope * i * dx for i in xrange(imax)]
        slope = umax / (imax * dx - self.length)
        q0 += [slope * (i * dx - self.length)
               for i in xrange(imax, self.n_modes + 2)]
        return npw.asrealarray(q0[1:-1])

    def compute_initial_state_modal(self, imax, umax):
        """Set initial positions of the string,
        assuming a triangular shape, with u[imax] = umax
        and modal form.
        """
        q0 = self._compute_initial_state_std(imax, umax)
        q0[...] = np.dot(self.s_mat.T, q0)
        coeff = self.length / (self.n_modes + 1)
        q0 *= coeff
        return npw.asrealarray(q0)

    def eigenfreq(self, j):
        """Compute eigenfrequency number j
        """
        return 0.5 * j * self.c0 / self.length * \
            math.sqrt(1 + self.stiffness_coeff * j)

    def compute_linear_coeff(self):
        """Compute M, K and C operators of the Lagrangian DS
        K = Omega^2
        C = 2.Gamma
        M = I
        """
        mass = np.ones(self.n_modes, dtype=np.float64)

        # --- Omega^2 vector ---
        omega2 = npw.asrealarray([1. + self.stiffness_coeff * (j ** 2)
                                  for j in xrange(1, self.n_modes + 1)])
        indices = np.arange(1, self.n_modes + 1)
        coeff = (indices * math.pi * self.c0 / self.length) ** 2
        omega2 *= coeff
        stiffness_mat = npw.asrealarray(omega2)

        # 2.S.Gamma.S-1
        sigma = self.compute_damping(np.sqrt(omega2) / (2. * math.pi))
        damping_mat = npw.asrealarray(sigma)
        damping_mat *= 2.
        return mass, stiffness_mat, damping_mat

    def compute_s_mat(self):
        """Compute 'S' matrix (normal modes)
        """
        indices = np.arange(1, self.n_modes + 1)
        row = indices.reshape(1, self.n_modes)
        s_mat = npw.asrealarray(row * row.T)
        s_mat *= (math.pi / self.n_modes)
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


class Fret(sk.Interaction):
    """Build a fret as a siconos interaction:

    * Newton impact nonsmooth law
    * Lagrangian Linear Time-invariant relation
    """
    def __init__(self, string, position=None, restitution_coeff=1.):
        """
        Parameters
        ----------
        string : :class:`StringDS`
            the dynamical system linked to this interaction
        positon : tuple or list
            position[0] = horizontal index position of the contact point
            position[1] = vertical position of the fret
        restitution_coeff : double, optional
            coefficient of restitution, default=1.
        """
        # siconos relation
        self.position = position[1] * npw.ones(1)
        self.contact_index = position[0]
        hmat = npw.zeros((1, string.n_modes))
        hmat[...] = string.s_mat[self.contact_index, :]
        e = restitution_coeff
        nslaw = sk.NewtonImpactNSL(e)
        relation = sk.LagrangianLinearTIR(hmat, -self.position)
        super(Fret, self).__init__(nslaw, relation)
