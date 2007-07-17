# Siconos-sample version 2.1.1, Copyright INRIA 2005-2007.
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.	
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# Siconos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY vincent.acary@inrialpes.fr 
#	
# SICONOS - PYTHON : bouncing ball example

import pySiconos;
from numpy import *
import Gnuplot
gnuplot = Gnuplot.Gnuplot(debug=1);
gnuplot.title("Python BouncingBall sample");

model = pySiconos.Model("./BallTS.xml");
print "*** BallTS.xml loaded ***";
s = pySiconos.TimeStepping.convert(model.getSimulationPtr());	
print "simulation initialization ..."; 
s.initialize();
print " **** the simulation is ready ****"; 
t = s.getTimeDiscretisationPtr();
k = 0;
N = t.getNSteps();

# DS
ball = pySiconos.LagrangianLinearTIDS.convert(model.getNonSmoothDynamicalSystemPtr().getDynamicalSystemPtr(0));

m = pySiconos.SimpleMatrix(N+1, 4);
m.zero()
#time
tmp = k * t.getH();
m.setValue(k, 0, tmp);
# position
val = (ball.getQPtr()).getValue(0);

m.setValue(k, 1, val);
# position
val = (ball.getVelocityPtr()).getValue(0);
m.setValue(k, 2, val);
val = (model.getNonSmoothDynamicalSystemPtr().getInteractionPtr(0).getLambdaPtr(1)).getValue(0);
m.setValue(k, 3, val);
# _________________________________________________

temps = pySiconos.doubleVector(N+1);
position = pySiconos.doubleVector(N+1);
plot = 0

while k < N :
	k = k+1;
	s.computeOneStep();
#	// Trace Values
#	//time
	m.setValue(k, 0, k*t.getH());
	temps[k] = k*t.getH();
	
#	// position
	m.setValue(k, 1, (ball.getQPtr()).getValue(0));
	position[k] = (ball.getQPtr()).getValue(0);
	
#	// position
	m.setValue(k, 2, (ball.getVelocityPtr()).getValue(0));						
	m.setValue(k, 3, (model.getNonSmoothDynamicalSystemPtr().getInteractionPtr(0).getLambdaPtr(1)).getValue(0));    
	if plot == 1999 :
		d1=Gnuplot.Data(temps, position, title="Ball position");
		gnuplot.plot(d1);
		plot = 0;
	else :
		plot = plot + 1
	s.nextStep();
print 'iterations  done: ', k;
io = pySiconos.ioMatrix("result.dat", "ascii");
io.write(m,"noDim");
#model.saveToXMLFile("./BouncingBall_TIDS.xml.output");

print '## END ##';
