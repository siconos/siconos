"""Fret (obstacle) definition as
a siconos interaction
"""

import siconos.kernel as sk
import numpywrappers as npw
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from scipy import signal


class Fret(sk.Interaction):
    """Build a fret as a siconos interaction
    """
    def __init__(self, string, position=None, restitution_coeff=1.):
        """"
        """
        # siconos relation
        self.position = position[1] * npw.ones(1)
        self.contact_index = position[0]
        if string.modal_form:
            coeff = 1.  # string.length / string._N
            if string.use_sparse:
                hmat = sk.SimpleMatrix(1, string.ndof, sk.SPARSE, string.ndof)
                for i in xrange(string.ndof):
                    val = coeff * string.s_mat[self.contact_index, i]
                    hmat.setValue(0, i, val)
            else:
                hmat = npw.zeros((1, string.ndof))
                hmat[...] = string.s_mat[self.contact_index, :]
                hmat *= coeff
        else:
            hmat = npw.zeros((1, string.ndof))
            hmat[0, self.contact_index] = 1.

        e = restitution_coeff
        nslaw = sk.NewtonImpactNSL(e)
        relation = sk.LagrangianLinearTIR(hmat, -self.position)
        super(Fret, self).__init__(nslaw, relation)


class Guitar(sk.Model):
    """DS (strings) and interaction (frets)
    'assembly' to build a NSDS
    """

    def __init__(self, strings, frets, time_range, fs):
        """
        """
        super(Guitar, self).__init__(time_range[0], time_range[1])
        if not isinstance(strings, list):
            strings = [strings]
        if not isinstance(frets, list):
            frets = [frets]
        self.modal_form = strings[0].modal_form
        self.frets = frets
        self.strings = strings
        for string in strings:
            self.nonSmoothDynamicalSystem().insertDynamicalSystem(string)
            # link the interaction(s) and the dynamical system(s)
            for fret in frets:
                self.nonSmoothDynamicalSystem().link(fret, string)

        # -- Simulation --
        # (1) OneStepIntegrators
        theta = 0.5001
        osi = sk.MoreauJeanOSI(theta)
        t0 = time_range[0]
        tend = time_range[1]
        # (2) Time discretisation --
        self.fs = fs  # 1960
        self.time_step = 1. / fs
        t = sk.TimeDiscretisation(t0, self.time_step)
        self.nb_time_steps = (int)((tend - t0) / self.time_step)
        # (3) one step non smooth problem
        osnspb = sk.LCP()
        # (4) Simulation setup with (1) (2) (3)
        self.simu = sk.TimeStepping(t, osi, osnspb)

        # simulation initialization
        self.setSimulation(self.simu)
        self.initialize()
        ndof = strings[0].ndof
        self.fret_position = frets[0].contact_index
        self.data = npw.zeros((self.nb_time_steps + 1, 3 + ndof))

    def save_state(self, k):
        """
        """
        self.data[k, 0] = self.simu.nextTime()
        if self.modal_form:
            self.data[k, 3:] = np.dot(self.strings[0].s_mat,
                                      self.strings[0].q())
        else:
            self.data[k, 3:] = self.strings[0].q()
        self.data[k, 1] = self.frets[0].y(0)
        self.data[k, 2] = self.frets[0].lambda_(1)

    def plot_state(self, nfig=1, pdffile=None):
        """

        Parameters
        ----------
        nfig : int
            figure number
        pdffile : string, optional
            output file name, if needed. Default=None
        """

        time = self.data[:, 0]
        dist = self.data[:, 1]
        lam = self.data[:, 2]
        qmax = self.data[:, 3 + self.fret_position]
        ndof = self.data.shape[1] - 3
        x = np.arange(ndof)
        plt.figure(nfig)
        plt.subplot(341)
        plt.plot(time, qmax)
        plt.title('displacements_' + str(self.modal_form))
        plt.axhline(self.frets[0].position, color='b', linewidth=3)
        plt.subplot(342)
        f, t, Sxx = signal.spectrogram(qmax, self.fs)
        plt.pcolormesh(t, f, Sxx)
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.title('dsp')
        plt.subplot(343)
        plt.plot(time, dist)
        plt.axhline(0, color='b', linewidth=3)
        plt.title('distance')
        plt.subplot(344)
        plt.plot(time, lam)
        plt.title('percussion')
        plt.subplot(345)
        plt.plot(x, self.data[0, 3:])
        plt.title('mode, t=0')
        plt.subplot(346)
        tint = self.nb_time_steps / 5
        plt.plot(x, self.data[tint, 3:])
        plt.title('mode, t1')
        plt.subplot(347)
        plt.plot(x, self.data[2 * tint, 3:])
        plt.title('mode, t2')
        plt.subplot(349)
        plt.plot(x, self.data[3 * tint, 3:])
        plt.title('mode, t3')
        plt.subplot(3, 4, 10)
        plt.plot(x, self.data[4 * tint, 3:])
        plt.title('mode, t4')
        plt.subplot(3, 4, 11)
        plt.plot(x, self.data[-1, 3:])
        plt.title('mode, t5')
        #if pdffile is not None:
        #    plt.savefig(pdffile, format='pdf')
        #plt.show()
        return plt

    def plot_modes(self, movie_name):
        """
        """

        fig = plt.figure()
        string = self.strings[0]
        ndof = string.ndof
        length = string.length
        ax = plt.axes(xlim=(0, length), ylim=(-0.01, 0.01))
        line, = ax.plot([], [], lw=2)

        # initialization function: plot the background of each frame
        def init():
            line.set_data([], [])
            return line,

        def animate(i):
            x = np.linspace(0., length, ndof)
            y = self.data[i, 3:]
            line.set_data(x, y)
            return line,

        #call the animator.
        # blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=self.nb_time_steps,
                                       interval=20, blit=True)
        anim.save(movie_name, fps=30, extra_args=['-vcodec', 'libx264'])
