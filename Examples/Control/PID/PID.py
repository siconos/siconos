from Siconos.Kernel import *
from matplotlib.pyplot import subplot, title, plot, grid, show
from numpy import array, eye, empty, zeros, savetxt
t0 = 0.0   # start time
T = 10      # end time
T = 100
h = 0.05;
xFinal = 0
theta = 0.5
A = zeros((2,2))
A[0,1] = 1

B = zeros(2)
x0 = array([10.,10.])
doubleIntegrator = FirstOrderLinearTIDS(x0,A,B)
process = Model(t0, T)
process.nonSmoothDynamicalSystem().insertDynamicalSystem(doubleIntegrator)
OSI = Moreau(doubleIntegrator,theta)
t = TimeDiscretisation(t0, h)
tSensor = TimeDiscretisation(t0, h)
tActuator = TimeDiscretisation(t0, h)
s = TimeStepping(t, 0)
s.insertIntegrator(OSI)
control = ControlManager(process)

C = array([1., 0], ndmin=2)
D = array([0, 0], ndmin=2)
sens = linearSensor(100, tSensor, process, C,D)
K = [.25, .125, 2]

# this fails !
control.addSensorPtr(sens)
