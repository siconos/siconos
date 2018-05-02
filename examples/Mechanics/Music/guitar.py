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
import h5py
import scipy.io
import os

class StringDS(sk.LagrangianLinearDiagonalDS):
    """Formulation for string-like structures,
    Based on LagrangianLinearDiagonalDS from Siconos

    Operators formulas from JSV Issanchou paper.
    """
    __damping_parameters_names = ['eta_air', 'rho_air', 'delta_ve', '1/qte']
    CLAMPED = 1

    def __init__(self, ndof, geometry_and_material,
                 max_coords=None, damping_parameters=None,
                 matlab_input=None):
        """Build string

        Parameters
        ----------
        ndof : int
            system size
        geometry_and_material : dictionnary
            description of the string, keys must be:
            ['length', 'diameter', 'B', 'tension', 'density']
        max_coords : tuple of double
            (umax, vmax), coordinates of the string maximum position
            at initial time (assuming triangular shape).
            vmax: along the neck, umax, vertical, displacement
        damping_parameters: dictionnary, optional
            constant parameters related to damping.
            keys must be ['eta_air', 'rho_air', 'delta_ve', '1/qte']
        matlab_input: string
            matlab file radix, used to read freq and damping values
            read from radix_frequs.mat and radix_amortissements.mat

        Notes
        -----
        * DS.stiffness = diag(omega**2)
        * DS.damping  = diag(2 * sigma)
        * DS.mass = identity

        """
        self.n_modes = ndof
        # B is a function of Young Modulus (E),
        # moment of inertia ...
        # B = pi^2EI / (TL^2)
        self.length = geometry_and_material['length']
        self.space_step = self.length / (self.n_modes + 1)
        self.x = np.linspace(0, self.length, ndof + 2)
        self.x = self.x[1:-1]
        self.s_mat = self.compute_s_mat()


        # - Compute operators (K, C, M)

        assert damping_parameters is not None or matlab_input is not None
        
        # -- First case : read material/geom parameters and compute operators of the DS --
        if damping_parameters is not None:
            assert matlab_input is None
            self.damping_parameters = damping_parameters
            msg = 'StringDS : missing parameter value for damping.'
            for name in self.__damping_parameters_names:
                assert name in self.damping_parameters.keys(), msg

            self.c0 = math.sqrt(self.tension / self.density)
            self.diameter = geometry_and_material['diameter']
            self.density = geometry_and_material['density']
            self.stiffness_coeff = geometry_and_material['B']
            self.tension = geometry_and_material['tension']

            stiffness_mat, damping_mat = self.compute_linear_coeff()
        
        # -- Second case : read material/geom parameters and compute operators of the DS --
        else:
            assert matlab_input is not None
            current_path = os.path.dirname(os.path.realpath(__file__))
            matlab_input = os.path.join(current_path, matlab_input)
            stiffness_mat, damping_mat = self.read_linear_coeff(matlab_input)
        self.max_coords = max_coords
        self.matlab_input = matlab_input
        # - Build siconos dynamical system -
        # -- initial conditions --
        qmod0 = self.compute_initial_state_modal(max_coords)
        v0 = npw.zeros_like(qmod0)
        super(StringDS, self).__init__(qmod0, v0, stiffness_mat, damping_mat)
        self.ic_status = self.check_initial_conditions()

    def check_initial_conditions(self):
        """Compare ds initial conditions with those
           read from matlab file.
        """
        # index and value of maximum of u0
        (val, ind) = (self.u0.max(), self.u0.argmax())
        assert np.allclose(self.x[ind], self.max_coords[1])
        assert np.allclose(val, self.max_coords[0])
        
        ic_filename = self.matlab_input + '_q2.mat'
        if os.path.exists(ic_filename):
            q0 = scipy.io.loadmat(ic_filename)['q2'][:, 0]
            return np.allclose(q0, self.q0(), atol=1e-6)
        else:
            return True
        
    def _compute_initial_state_std(self, max_coords):
        """Set initial positions of the string,
        assuming a triangular shape, with u[imax] = umax.
        """
        umax = max_coords[0]
        # index of the max
        dx = self.space_step
        imax = int(round(max_coords[1] / dx))
        slope = umax / (imax * dx)
        u0 = [slope * i * dx for i in range(imax)]
        slope = umax / (imax * dx - self.length)
        u0 += [slope * (i * dx - self.length)
               for i in range(imax, self.n_modes + 2)]
        return npw.asrealarray(u0[1:-1])

    def compute_initial_state_modal(self, max_coords):
        """Set initial positions of the string,
        assuming a triangular shape, with u[imax] = umax
        and modal form.
        """
        if max_coords is None:
            assert self.matlab_input is not None
            inputfile = self.matlab_input + '_q2.mat'
            q0 = scipy.io.loadmat(inputfile)['q2'][:, 0].copy()
            self.u0 = np.dot(self.s_mat, q0)
            self.max_coords = (self.u0.max(), self.x[self.u0.argmax()])
            return q0
        else:
            self.u0 = self._compute_initial_state_std(max_coords)
            q0 = np.dot(self.s_mat.T, self.u0)
            coeff = self.length / (self.n_modes + 1)
            q0 *= coeff
            return npw.asrealarray(q0)

    def eigenfreq(self, j):
        """Compute eigenfrequency number j
        """
        return 0.5 * j * self.c0 / self.length * \
            np.sqrt(1 + self.stiffness_coeff * j)

    def read_linear_coeff(self, radix):
        """Read K and C operators of the Lagrangian DS
        K = Omega^2
        C = 2.Gamma
        from matlab files
        """
       
        freq_file = radix + '_frequs.mat'
        damp_file = radix + '_amortissements.mat'
        nuj = scipy.io.loadmat(freq_file)['frequs'][:, 0]
        assert nuj.size == self.n_modes
        omega2 = (2. * np.pi * nuj) ** 2
        stiffness_mat = npw.asrealarray(omega2)
        if os.path.exists(damp_file):
            sigmas = scipy.io.loadmat(damp_file)['sig0'][:, 0]
            assert sigmas.size == self.n_modes
            damping_mat = 2. * npw.asrealarray(sigmas)
        else:
            print("Warning: no input for damping. Set damping matrix to zero.")
            damping_mat = None
        return stiffness_mat, damping_mat

    def compute_linear_coeff(self):
        """Compute M, K and C operators of the Lagrangian DS
        K = Omega^2
        C = 2.Gamma
        M = I
        """
        #mass = np.ones(self.n_modes, dtype=np.float64)
        # --- Omega^2 vector ---
        # Ref : Clara's phd manuscript, eq 2.4
        omega2 = npw.asrealarray([1. + self.stiffness_coeff * (j ** 2)
                                  for j in range(1, self.n_modes + 1)])
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
        # Compute quality factor as defined in eq 2.20 from Clara's manuscript.
        #
        eta_air = self.damping_parameters['eta_air']
        rho_air = self.damping_parameters['rho_air']
        delta_ve = self.damping_parameters['delta_ve']
        r = np.sqrt(math.pi * eta_air * rho_air * nu)
        r *= self.diameter
        r += eta_air
        r *= 2. * math.pi
        indices = np.arange(1, self.n_modes + 1)
        q_air_inv = r / (2. * math.pi * self.density * nu)
        coeff = 0.5 * indices * self.c0 / self.length
        q_air_inv *= coeff
        q_air_inv /= nu
        q_ve_inv = 4 * self.length ** 2 * self.stiffness_coeff * self.density
        q_ve_inv *= delta_ve
        nu_0 = 0.5 * indices * self.c0 / self.length
        q_ve_inv *= nu_0 ** 3
        q_ve_inv /= nu
        q_ve_inv /= self.tension
        quality_factor_inv = q_ve_inv + q_air_inv
        quality_factor_inv += self.damping_parameters['1/qte']
        return quality_factor_inv * math.pi * nu


class Fret(sk.Interaction):
    """Build a fret as a siconos interaction:

    * Newton impact nonsmooth law
    * Lagrangian Linear Time-invariant relation
    """
    def __init__(self, string, contact_positions, restitution_coeff):
        """
        Parameters
        ----------
        string : :class:`StringDS`
            the dynamical system linked to this interaction
        contact_positons : tuple
            contact_positions[0] = horizontal index position
             of the contact point
            contact_positions[1] = vertical position of the fret
        restitution_coeff : double
            coefficient of restitution
        """
        # vertical positions of the contact points
        hmat = npw.zeros((1, string.n_modes))
        dx = string.space_step
        hmat[0, :] = string.s_mat[contact_positions[0], :]
        # compute contact point horizontal position
        self.contact_pos = dx * contact_positions[0]
        # set contact index (mind that boundary points
        # are not included in fret/ds)
        self.contact_index = contact_positions[0]
        # Build nslaw, relation and interaction
        e = restitution_coeff
        nslaw = sk.NewtonImpactNSL(e)
        dist = -contact_positions[1]
        relation = sk.LagrangianLinearTIR(hmat, [dist])
        super(Fret, self).__init__(nslaw, relation)


class Guitar(object):
    """DS (strings) and interaction (frets)
    'assembly' to build a NSDS
    """

    # Available settings for interactions outputs
    __authorized_outputs__ = [1, 2, 3]
    
    def __init__(self, strings_and_frets, time_range, fs, output_freq=1,
                 interactions_output=0,
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
        output_freq: int, optional
            output frequency for times steps
        interactions_output: int, optional
            0: save nothing, 1: only y, 2: y + lambda(vel) 3: y + ydot + lambda
            default = 0
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
        # iterate through ds and interactions in the input dictionnary
        # - insert ds into the nonsmooth dynamical system
        # - link each ds to its interaction
        self.nsds = sk.NonSmoothDynamicalSystem(time_range[0], time_range[1])
        if None in strings_and_frets:
            # No interactions in the model
            ds = strings_and_frets[None]
            assert isinstance(ds, StringDS)
            self.nsds.insertDynamicalSystem(ds)
        else:
            self.nsds.insertDynamicalSystem(list(strings_and_frets.values())[0])
            for interaction in strings_and_frets:
                ds = strings_and_frets[interaction]
                assert isinstance(interaction, Fret)
                assert isinstance(ds, StringDS)
                #self.nsds.insertDynamicalSystem(ds)
                # link the interaction and the dynamical system
                self.nsds.link(interaction, ds)
        self.strings_and_frets = strings_and_frets
        # -- Simulation --
        moreau_bilbao = sk.MoreauJeanBilbaoOSI()
        #moreau_jean = sk.MoreauJeanOSI(0.500001)
        #default_integrator = moreau_bilbao
        #if default_integrator == 'MoreauJean':
        #    default_integrator = moreau_jean
        
        # (2) Time discretisation --
        t0 = time_range[0]
        tend = time_range[1]
        self.fs = fs
        self.time_step = 1. / fs
        t = sk.TimeDiscretisation(t0, self.time_step)
        self.nb_time_steps = (int)((tend - t0) * fs)
        self.output_freq = output_freq
        self.nb_time_steps_output = (int)(self.nb_time_steps / output_freq)
        # (3) one step non smooth problem
        self.osnspb = sk.LCP()
        # (4) Simulation setup with (1) (2) (3)
        self.default_integrator = moreau_bilbao #default_integrator
        #self.simulation = sk.TimeStepping(self.nsds, t, default_integrator, self.osnspb)
        self.simulation = sk.TimeStepping(self.nsds, t, moreau_bilbao, self.osnspb)
        # if integrators is not None:
        #     for ds in integrators:
        #         self.simulation.prepareIntegratorForDS(integrators[ds], ds, self, t0)

        # internal buffers, used to save data for each time step.
        # A dict of buffers to save ds variables for all time steps
        self.data_ds = {}
        for ds in self.strings_and_frets.values():
            ndof = ds.dimension()
            self.data_ds[ds] = npw.zeros((ndof, self.nb_time_steps_output + 1))
        if None in self.strings_and_frets:
            self.strings_and_frets = {}

        # A dict of buffers to save interactions variables for all time steps
        self.data_interactions = {}
        self.save_interactions = interactions_output > 0
        self.interactions_output = interactions_output
        if self.save_interactions:
            for interaction in self.strings_and_frets:
                self.data_interactions[interaction] = []
                for i in range(interactions_output):
                    self.data_interactions[interaction].append(
                        npw.zeros(self.nb_time_steps_output + 1))
        # time instants
        self.time = npw.zeros(self.nb_time_steps_output + 1)

        # saved data state (modal or nodal)
        self.modal_values = True
        
    def save_ds_state_modal(self, k, ds):
        """Save ds modal positions, at iteration k

        Parameters
        ----------
        k : int
            current iteration number
        ds : StringDS
            dynamical system of interest
        """
        self.data_ds[ds][:, k] = ds.q()

    def convert_modal_output(self, ds):
        """Post-processing.
        Recover nodal values from modal outputs.

        Parameters
        ----------
        ds : StringDS
            dynamical system of interest
        """
        if self.modal_values:
            self.data_ds[ds][:, :] = np.dot(ds.s_mat, self.data_ds[ds])
            self.modal_values = False


    def plot_traj(self, ds, dof, filename=None, iplot=0, ground=None):
        """Plot collected data (positions ...) of a dynamical system

        Parameters
        ----------
        ds : StringDS
            dynamical system of interest
        dof : int
            index of dof to be plotted.
        filename : string, optional
            name of the output file
        iplot : int
            parent figure number
        """
        self.convert_modal_output(ds)
        data = self.data_ds[ds]
        # current interaction
        # number of contact points
        # distance(s) string/fret
        ndof = ds.dimension()
        x = np.linspace(0, ds.length, ndof + 2)
        x = x[1:-1]
        # Plot string displacements, at contact points, according to time
        plt.figure(iplot, figsize=(17, 8))
        leg = []
        plt.subplot(2, 2, 1)
        plt.plot(self.time, data[dof, :])
        if ground is not None:
            plt.plot((self.time[0], self.time[-1]),(ground, ground), '-')
        plt.subplot(2, 2, 2)
        plt.plot(self.time, data[dof, :], 'x-')
        plt.xlim(0, 0.008)
        #plt.xlim(0.5, 1.)
        if ground is not None:
            plt.plot((self.time[0], self.time[-1]),(ground, ground), '-x')

        plt.subplot(2, 2, 3)
        plt.plot(self.time, data[dof, :])
        plt.xlim(0.05, 0.07)
        #plt.xlim(1.45, 1.55)
        if ground is not None:
            plt.plot((self.time[0], self.time[-1]),(ground, ground), '-x')

        plt.subplot(2, 2, 4)
        plt.plot(self.time, data[dof, :])
        plt.xlim(0.35, 0.377)
        #plt.xlim(2.5, 3.5)
        if ground is not None:
            plt.plot((self.time[0], self.time[-1]),(ground, ground), '-x')
        leg.append('x = ' + str(x[dof]))
        plt.legend(leg)
        plt.suptitle('displacements = f(time) at x='+str(x[dof]))
        if filename is not None:
            plt.savefig(filename)
        return plt

    def plot_modes(self, ds, times=None,
                   plot_shape=None, filename=None, iplot=1):
        """Plot collected data (positions ...) of a dynamical system

        Parameters
        ----------
        ds : StringDS
            dynamical system of interest
        plot_shape : tuple
            subplot (i.e. grid of figures) shape
        filename : string, optional
            name of the output file
        iplot : int
            parent figure number
        """
        if times is None:
            if plot_shape is None:
                plot_shape = [2, 4]

            # Split time range using the number of required figures
            nb_points = plot_shape[0] * plot_shape[1]
            plot_x = plot_shape[0]
            plot_y = plot_shape[1]
            nb_points = plot_x * plot_y
            freq = self.nb_time_steps_output // nb_points
            time_ind = np.arange(0, self.nb_time_steps_output, freq)
            time_ind[-1] = -1

        else:
            time_ind = times
            nb_points = len(times)
            if plot_shape is None:
                plot_shape = (nb_points // 2 + nb_points % 2, 2)
            plot_x = plot_shape[0]
            plot_y = plot_shape[1]

        ndof = ds.dimension()
        plt.figure(iplot, figsize=(17, 8))
        # plt.subplot(342)
        # #f, t, Sxx = signal.spectrogram(pos[0], self.fs)
        # #plt.pcolormesh(t, f, Sxx)
        # plt.ylabel('Frequency [Hz]')
        # plt.xlabel('Time [sec]')
        # plt.title('dsp')
        # plt.subplot(343)
        # output frequency, for modes
        
        self.convert_modal_output(ds)
        data = self.data_ds[ds]
        interactions = self.interactions_linked_to_ds(ds)
        nbc = len(interactions)
        pos = 1
        ymin = data[:, time_ind].min()
        for k in range(nb_points):
            plt.subplot(plot_x, plot_y, pos)
            plt.plot(ds.x, data[:, time_ind[k]])
            plt.title('mode, t=' + str(self.time[time_ind[k]]))
            ylimits = (ymin - 0.2 * abs(ymin),
                       1.1 * data[:, time_ind[k]].max())
            for ic in range(nbc):
                vpos = -interactions[ic].relation().e()[0]
                plt.plot((interactions[ic].contact_pos,
                          interactions[ic].contact_pos),
                         (2. * vpos, vpos),
                         'o-')
            plt.ylim(ylimits)
            pos += 1
        plt.subplots_adjust(hspace=0.8)
        if filename is not None:
            plt.savefig(filename)

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
        assert self.save_interactions, 'Interactions output is not enabled.'
        data = self.data_interactions[interaction]
        # current interaction
        # number of contact points
        nbc = interaction.dimension()
        # distance(s) string/fret
        dist = data[:, :nbc]
        vel = data[:, nbc:2 * nbc]
        # reaction(s) (impulse) at contact(s)
        lam = data[:, 2 * nbc:]
        plt.figure(nfig, figsize=(17, 8))
        plt.subplot(131)
        plt.plot(self.time, dist)
        plt.axhline(0, color='b', linewidth=3)
        plt.title('distance')
        plt.subplot(132)
        plt.plot(self.time, vel)
        plt.title('velo')
        plt.subplot(133)
        plt.plot(self.time, lam)
        plt.title('percussion')
        return plt

    def make_movie(self, ds, movie_name, sampling=100):
        """Create animation from simulation results,
        for a given ds.
        """
        self.convert_modal_output(ds)
        data = self.data_ds[ds]
        ylimits = (data.min() - 0.2 * abs(data.min()),
                   1.1 * data.max())

        interactions = self.interactions_linked_to_ds(ds)
        nbc = len(interactions)
        fig = plt.figure()
        ndof = ds.dimension()
        length = ds.length
        ax = plt.axes(xlim=(0, length), ylim=ylimits)
        line, = ax.plot([], [], lw=2)
        for ic in range(nbc):
            vpos = -interactions[ic].relation().e()[0]
            plt.plot((interactions[ic].contact_pos,
                      interactions[ic].contact_pos),
                     (2. * vpos, vpos),
                     'o-')

        # initialization function: plot the background of each frame
        def init():
            line.set_data([], [])
            return line,

        def animate(i):
            y = self.data_ds[ds][:, i]
            x = ds.x
            line.set_data(x, y)
            return line,

        #call the animator.
        # blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(
            fig, animate, init_func=init,
            frames=range(0,self.nb_time_steps_output,sampling),
            blit=True)
        anim.save(movie_name, fps=100, extra_args=['-vcodec', 'libx264'])

    def contactogram(self, ds, nfig=12):
        """Plot contact times on each fret
           for a given string

        Parameters
        ----------
        ds : StringDS
            dynamical system of interest

        """
        assert self.save_interactions, 'Interactions output is not enabled.'
        interactions = self.interactions_linked_to_ds(ds)
        nb_inter = len(interactions)
        print('nb contacts : ', nb_inter)
        plt.figure(nfig, figsize=(17, 8))
        for ic in range(nb_inter):
            inter = interactions[ic]
            #nbc = inter.dimension()
            # find lambda > 0 to identify contact times
            contact_indices = np.where(
                self.data_interactions[inter][1][:] > 0.)# 1e-8)
            nbcontacts = len(contact_indices[0])
            pos = inter.contact_pos
            plt.subplot(1,2,1)
            plt.plot(self.time[contact_indices[0]], [pos, ] * nbcontacts, 'o')
            plt.subplot(1,2,2)
            plt.plot(self.time[contact_indices[0]], [pos, ] * nbcontacts, 'o')
            plt.xlim(0.0038, 0.00425 )
            # plt.subplot(3,2,1)
            # plt.plot(self.time[contact_indices[0]], [pos, ] * nbcontacts, 'o')
            # plt.subplot(3,2,2)
            # plt.plot(self.time[contact_indices[0]], [pos, ] * nbcontacts, 'o')
            # plt.xlim(0.00367, 0.0040 )
            # plt.subplot(3,2,3)
            # plt.plot(self.time[contact_indices[0]], [pos, ] * nbcontacts, 'o')
            # plt.xlim(0.0115, 0.0125 )
            # #plt.ylim(0.18, 0.33 )
            # plt.subplot(3,2,4)
            # plt.plot(self.time[contact_indices[0]], [pos, ] * nbcontacts, 'o')
            # plt.xlim(0.02125, 0.0215 )
            # plt.ylim(0.1, 0.21 )
            # plt.subplot(3,2,5)
            # plt.plot(self.time[contact_indices[0]], [pos, ] * nbcontacts, 'o')
            # plt.xlim(0.03125, 0.0315 )
            # plt.ylim(0.12, 0.2 )
        #plt.yticks(np.arange(0, nb_inter, ))
        #    plt.xlim(0, self.time[-1])
        plt.xlabel('time')
        plt.ylabel('frets positions')
