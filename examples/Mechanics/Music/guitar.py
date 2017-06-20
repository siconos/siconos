"""Definition of string dynamics (class StringDS)
and fret interaction (class Fret) as objects derived
from Siconos classes.
"""
import math
import numpy as np
import siconos.kernel as sk
# import scipy.sparse as scs
import numpywrappers as npw
import matplotlib.pyplot as plt
#from scipy import signal
from matplotlib import animation


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
        stiffness_mat, damping_mat = self.compute_linear_coeff()
        q0 = self.compute_initial_state_modal(imax, umax)
        v0 = npw.zeros_like(q0)
        super(StringDS, self).__init__(q0, v0, stiffness_mat, damping_mat)

    def _compute_initial_state_std(self, imax, umax):
        """Set initial positions of the string,
        assuming a triangular shape, with u[imax] = umax.
        """
        dx = self.length / (self.n_modes + 1)
        slope = umax / (imax * dx)
        u0 = [slope * i * dx for i in xrange(imax)]
        slope = umax / (imax * dx - self.length)
        u0 += [slope * (i * dx - self.length)
               for i in xrange(imax, self.n_modes + 2)]
        return npw.asrealarray(u0[1:-1])

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
        #mass = np.ones(self.n_modes, dtype=np.float64)
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
        return stiffness_mat, damping_mat

    def compute_s_mat(self):
        """Compute 'S' matrix (normal modes)
        """
        nn = self.n_modes + 1
        indices = np.arange(1, nn)
        row = indices.reshape(1, self.n_modes)
        s_mat = npw.asrealarray(row * row.T)
        s_mat *= (math.pi / nn)
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
    def __init__(self, string, contact_positions, restitution_coeff=1.):
        """
        Parameters
        ----------
        string : :class:`StringDS`
            the dynamical system linked to this interaction
        contact_positons : tuple
            contact_positions[0] = horizontal index position
             of the contact point
            contact_positions[1] = vertical position of the fret
        restitution_coeff : double, optional
            coefficient of restitution, default=1.
        """
        # vertical positions of the contact points
        hmat = npw.zeros((1, string.n_modes))
        dx = string.length / (string.n_modes + 1)
        hmat[0, :] = string.s_mat[contact_positions[0] - 1, :]
        # compute contact point horizontal position
        self.contact_pos = dx * contact_positions[0]
        # set contact index (mind that boundary points
        # are not included in fret/ds)
        self.contact_index = contact_positions[0] - 1
        # Build nslaw, relation and interaction
        e = restitution_coeff
        nslaw = sk.NewtonImpactNSL(e)
        relation = sk.LagrangianLinearTIR(hmat, [-contact_positions[1]])
        super(Fret, self).__init__(nslaw, relation)


class Guitar(sk.Model):
    """DS (strings) and interaction (frets)
    'assembly' to build a NSDS
    """

    def __init__(self, strings_and_frets, time_range, fs,
                 integrators=None, default_integrator=None):
        """
        Parameters
        ----------
        strings_and_frets : python dictionnary
            keys = interactions (Fret). See notes below.
            values = dynamical systems (StringDS)
        time_range : list
            time range for the simulation.
        fs : double
            sampling frequency
        integrators : dictionnary, optional
            integrators[some_ds] = osi class name
            (MoreauJeanOSI or MoreauJeanBilbaoOSI).
            Default: Bilbao for all ds.
        default_integrator: sk.OneStepIntegrator
            default osi type (if integrators is not set
            or for ds not present in integrators).
            Default = Bilbao.

        Notes
        -----
        * strings_and_frets.keys() : interactions
        * strings_and_frets.values() : dynamical systems
        * so, strings_and_frets[some_interaction] returns the ds
        linked to some_interaction.
        * For integrators:
          - to use Bilbao for all ds : do not set
            integrators and default_integrator params
          - to choose the integrator: set default_integrator=OSI
            OSI being some siconos OneStepIntegrator class.
          - to explicitely set the integrator for each ds:
            set integrators dictionnary, with
            integrators[some_ds] = OSI (OneStepIntegrator class)
            all ds not set in integrators will be integrated
            with default_integrator.
        """
        # Build siconos model object (input = time range)
        super(Guitar, self).__init__(time_range[0], time_range[1])
        # iterate through ds and interactions in the input dictionnary
        # - insert ds into the nonsmooth dynamical system
        # - link each ds to its interaction
        nsds = self.nonSmoothDynamicalSystem()
        for interaction in strings_and_frets:
            ds = strings_and_frets[interaction]
            assert isinstance(interaction, Fret)
            assert isinstance(ds, StringDS)
            nsds.insertDynamicalSystem(ds)
            # link the interaction and the dynamical system
            nsds.link(interaction, ds)
        self.strings_and_frets = strings_and_frets
        # -- Simulation --
        moreau_bilbao = sk.MoreauJeanBilbaoOSI()
        moreau_jean = sk.MoreauJeanOSI(0.500001)
        default_integrator = moreau_bilbao
        if default_integrator is 'MoreauJean':
            default_integrator = moreau_jean

        # (2) Time discretisation --
        t0 = time_range[0]
        tend = time_range[1]
        self.fs = fs
        self.time_step = 1. / fs
        t = sk.TimeDiscretisation(t0, self.time_step)
        self.nb_time_steps = (int)((tend - t0) / self.time_step)
        # (3) one step non smooth problem
        osnspb = sk.LCP()
        # (4) Simulation setup with (1) (2) (3)
        self.simu = sk.TimeStepping(t, default_integrator, osnspb)
        if integrators is not None:
            self.simu.insertIntegrator(integrators[ds])
            for ds in integrators:
                nsds.setOSI(ds, integrators[ds])
        # simulation initialization
        self.setSimulation(self.simu)
        self.initialize()
        # internal buffers, used to save data for each time step.
        # A dict of buffers to save ds variables for all time steps
        self.data_ds = {}
        for ds in self.strings_and_frets.values():
            ndof = ds.dimension()
            self.data_ds[ds] = npw.zeros((self.nb_time_steps + 1, ndof))
        # A dict of buffers to save interactions variables for all time steps
        self.data_interactions = {}
        for interaction in self.strings_and_frets:
            nbc = interaction.getSizeOfY()
            self.data_interactions[interaction] = \
                npw.zeros((self.nb_time_steps + 1, 2 * nbc))
        # time instants
        self.time = npw.zeros(self.nb_time_steps + 1)

    def save_interaction_state(self, k, interaction):
        """Save ds positions, velocity,
        and contact points local variables.

        Parameters
        ----------
        k : int
            current iteration number
        interaction : Fret
            interaction of interest
        """
        nbc = interaction.getSizeOfY()
        self.data_interactions[interaction][k, 1:1 + nbc] = interaction.y(0)
        self.data_interactions[interaction][k, 1 + nbc:] = \
            interaction.lambda_(1)

    def save_ds_state(self, k, ds):
        """Save ds positions, velocity,
        and contact points local variables.

        Parameters
        ----------
        k : int
            current iteration number
        ds : StringDS
            dynamical system of interest
        """
        self.data_ds[ds][k, :] = np.dot(ds.s_mat, ds.q())

    def plot_ds_state(self, ds, indices=None,
                      nfig=1, pdffile=None):
        """Plot collected data (positions ...) of a dynamical system

        Parameters
        ----------
        ds : StringDS
            dynamical system of interest
        indices : list of int, optional
            indices (dof) to be plotted. If None
            plot all dof.
        nfig : int
            figure number
        pdffile : string, optional
            output file name, if needed. Default=None
        """
        data = self.data_ds[ds]
        # current interaction
        # number of contact points
        # distance(s) string/fret
        ndof = ds.dimension()
        x = np.linspace(0, ds.length, ndof + 2)
        x = x[1:-1]
        plt.figure(nfig, figsize=(17, 8))
        # Plot string displacements, at contact points, according to time
        plt.subplot(341)
        if indices is None:
            indices = np.arange(ndof)
            for ind in indices:
                plt.plot(self.time, data[:, ind])
            plt.title('displacements')
        else:
            iplot = 1
            for ind in indices:
                plt.subplot(3, 4, iplot)
                if iplot < 5:
                    iplot += 1
                # plot ind - 1 because boundaries points
                # are not included in the ds
                plt.plot(self.time, data[:, ind - 1])
                plt.title('displacements at dof ' + str(ind))

        # plt.subplot(342)
        # #f, t, Sxx = signal.spectrogram(pos[0], self.fs)
        # #plt.pcolormesh(t, f, Sxx)
        # plt.ylabel('Frequency [Hz]')
        # plt.xlabel('Time [sec]')
        # plt.title('dsp')
        # plt.subplot(343)
        # output frequency, for modes
        nb_points = 8
        time_ind = np.arange(nb_points) * self.nb_time_steps / nb_points
        time_ind[-1] = -1
        interactions = self.interactions_linked_to_ds(ds)
        nbc = len(interactions)
        for k in range(nb_points):
            plt.subplot(3, 4, 5 + k)
            plt.plot(x, data[time_ind[k], :])
            plt.title('mode, t=' + str(self.time[time_ind[k]]))
            for ic in range(nbc):
                plt.plot((interactions[ic].contact_pos,
                          interactions[ic].contact_pos),
                         (-1.e-4, -interactions[ic].relation().e()[0]),
                         'o-')
        plt.subplots_adjust(hspace=0.8)
        return plt

    def interactions_linked_to_ds(self, ds):
        """Return a list of all interactions linked to a given
        ds.
        """
        return [k for k in self.strings_and_frets
                if self.strings_and_frets[k] in [ds]]

    def plot_interaction(self, interaction, nfig=1):
        """Plot collected data (positions ...) of an interaction

        Parameters
        ----------
        interaction : Fret
             interaction of interest
        nfig : int
            figure number
        """
        data = self.data_interactions[interaction]
        # current interaction
        # number of contact points
        nbc = interaction.getSizeOfY()
        # distance(s) string/fret
        dist = data[:, :nbc]
        # reaction(s) (impulse) at contact(s)
        lam = data[:, nbc:2 * nbc]
        plt.figure(nfig, figsize=(17, 8))
        plt.subplot(121)
        plt.plot(self.time, dist)
        plt.axhline(0, color='b', linewidth=3)
        plt.title('distance')
        plt.subplot(122)
        plt.plot(self.time, lam)
        plt.title('percussion')
        return plt

    def plot_modes(self, ds, movie_name):
        """Create animation from simulation results,
        for a given ds.
        """
        interactions = self.interactions_linked_to_ds(ds)
        nbc = len(interactions)
  
        fig = plt.figure()
        ndof = ds.dimension()
        length = ds.length
        ax = plt.axes(xlim=(0, length), ylim=(-2.e-3, 2e-3))
        line, = ax.plot([], [], lw=2)
        for ic in range(nbc):
            plt.plot((interactions[ic].contact_pos,
                      interactions[ic].contact_pos),
                     (-1.e-4,
                      -interactions[ic].relation().e()[0]),
                     'o-')

        # initialization function: plot the background of each frame
        def init():
            line.set_data([], [])
            return line,

        def animate(i):
            x = np.linspace(0., length, ndof)
            y = self.data_ds[ds][i, :]
            line.set_data(x, y)
            return line,

        #call the animator.
        # blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=int(self.nb_time_steps / 20.),
                                       interval=20, blit=True)
        anim.save(movie_name, fps=30, extra_args=['-vcodec', 'libx264'])
