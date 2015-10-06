#! /usr/bin/env python


from numpy import sin, cos, pi, array
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

L1 = 1.0
L2 = 1.0

# create a time array from 0..100 sampled at 0.1 second steps
dt = 0.005
t = np.arange(0.0, 20, dt)

f = np.loadtxt('DoublePendulumResult.dat')

x1 = f[::10, 5]
y1 = f[::10, 6]

x2 = f[::10, 7]
y2 = f[::10, 8]

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template%(i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(f[:])/10),
    interval=25, blit=True, init_func=init)

#ani.save('double_pendulum.mp4', fps=15, clear_temp=True)
plt.show()
