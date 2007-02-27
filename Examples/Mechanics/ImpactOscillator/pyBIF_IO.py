# Siconos-sample version 2.0.1, Copyright INRIA 2005-2006.
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
# SICONOS - PYTHON : impact oscillator example

import pySiconos;
#import funct;
import math;
import time;
import string;
from Numeric import *
import Gnuplot
gnuplot = Gnuplot.Gnuplot(debug=1);
gnuplot.title("Python Impact Oscillator")

model = pySiconos.Model("./IO.xml");

#===============================


# # Oscillator
q0 = pySiconos.SimpleVector(2);
q0.zero();
q0.setValue(0, .5)

v0 = pySiconos.SimpleVector(2);
v0.zero();
#v0.setValue(0, .1)


parmax=2.6;
parmin=2.6;

no_par=1;

hh=(parmax-parmin)/no_par;

omega=parmin;
trans_iter=2;
per_iter=2;
no_iter = 20;


xp = pySiconos.SimpleMatrix(trans_iter, 1);
yp = pySiconos.SimpleMatrix(trans_iter, 1);

s = model.getSimulationPtr();	
print "simulation will be initialized"; 
s.initialize();

print " **** the simulation is ready ****"; 
#t = s.getTimeDiscretisationPtr();
#k = t.getK();
#N = t.getNSteps();

ds1 = pySiconos.LagrangianLinearTIDS.convert(model.getNonSmoothDynamicalSystemPtr().getDynamicalSystemPtr(0));
ds2 = pySiconos.LagrangianLinearTIDS.convert(model.getNonSmoothDynamicalSystemPtr().getDynamicalSystemPtr(1));


ImpOsc = pySiconos.LagrangianLinearTIDS.convert(model.getNonSmoothDynamicalSystemPtr().getDynamicalSystemPtr(0));

Tk = pySiconos.SimpleVector(no_iter+1);
t = s.getTimeDiscretisationPtr();
t.setTkPtr(Tk);


for ii in range(no_par+1):
	
	Tend = 2*math.pi/omega;

	q0 = pySiconos.SimpleVector(2);
	q0.setValue(0,1);
	q0.setValue(1,omega);

	v0 = pySiconos.SimpleVector(2);
	v0.zero();

	print "Iteration aa:", q0.getValue(1)

	ds1.setQ(q0);
	print "Iteration bb "
	ds1.setVelocity(v0);
	print "Iteration cc "


	N = t.getNSteps();

	t.setT(Tend);
	t.setH(Tend/no_iter);
	k = 0;
	t.setNSteps(no_iter);

	for tt in range(no_iter+1):
		(t.getTkPtr()).setValue(tt,tt*(t.getH()));

	print "Iteration A ", ii
	t.display();
	print "Iteration B ", ii
      	print "Iteration C ", ii
	t.display();
	print "Iteration D ", ii


	
	# The for/while-loop
	kk = 0;


	print 'adra1 = ', (ImpOsc.getQPtr()).getValue(0),(ImpOsc.getVelocityPtr()).getValue(0);
	print 'adra2 = ', (ImpOsc.getQPtr()).getValue(1),(ImpOsc.getVelocityPtr()).getValue(1);

	while kk < trans_iter:
#		print "Iteration A ", kk

		tk = pySiconos.intVector(2);
		tk[0] = 0; tk[1] = 1;

		t = s.getTimeDiscretisationPtr();	
		t.setK(0);
		model.setCurrentT(0);
		k++;
		N = t.getNSteps();
#		print " **** time management OK ****";


		while k < N :
#			print "SOLVER ITERATION1", k

	#		print 'NextStep done';
			k++;
			s.computeOneStep();
			print '1andra = ', (ImpOsc.getQPtr()).getValue(0),(ImpOsc.getVelocityPtr()).getValue(0), ' time = ',(t.getTkPtr()).getValue(k);
			print '2andra = ', (ImpOsc.getQPtr()).getValue(1),(ImpOsc.getVelocityPtr()).getValue(1), ' time = ',(t.getTkPtr()).getValue(k);
			s.nextStep();

		print 'Final1', (ImpOsc.getQPtr()).getValue(0),(ImpOsc.getVelocityPtr()).getValue(0), ' time = ',t.getT();

#		ImpOsc = pySiconos.LagrangianLinearTIDS.convert(model.getNonSmoothDynamicalSystemPtr().getDynamicalSystemPtr(0));

		xp.setValue(kk,0,(ImpOsc.getQPtr()).getValue(0));
		yp.setValue(kk,0,(ImpOsc.getVelocityPtr()).getValue(0));

#		print "SOLVER ITERATION AA", k

		q0.setValue(0, (ImpOsc.getQPtr()).getValue(0));
		v0.setValue(0, (ImpOsc.getVelocityPtr()).getValue(0));
		ds1.setQ(q0);
		ds1.setVelocity(v0);		

		kk = kk+1;

#		print 'KK = ',kk
#		time.sleep(3)
		
	kkk=0;
	points=[];

	xp.write("XP0001.dat", "ascii");		
	yp.write("YP0001.dat", "ascii");		

	while kkk < per_iter:

		s = model.getSimulationPtr();
		#print "simulation will be initialized";
		#s.initialize();
		#print " **** the simulation is ready ****";

#		t = s.getTimeDiscretisationPtr();
		model.setCurrentT(0);
		k = 0;
		N = t.getNSteps();
#		print " **** time management OK ****";

		while k < N :
	#		print "SOLVER ITERATION1", k
	#		print 'NextStep done';
			k++;
			s.computeOneStep();
			s.nextStep();
		
		value=(ImpOsc.getQPtr()).getValue(0);

		jj=0;
		
		if len(points)==0:
			points=points+[value];
		else:
			while jj<len(points):
				if abs(points[jj]-value)>10**(-10):
					jj=len(points);
					points=points+[value];
				jj=jj+1;


		kkk = kkk+1;
		
		print 'Points:', points
		print 'time:', k*t.getH();

		print 'Final2', (ImpOsc.getQPtr()).getValue(0),(ImpOsc.getVelocityPtr()).getValue(0);
	#	time.sleep(3)

		q0.setValue(0, (ImpOsc.getQPtr()).getValue(0));
		v0.setValue(0, (ImpOsc.getVelocityPtr()).getValue(0));
		ds1.setQ(q0);
		ds1.setVelocity(v0);


	omega=omega+hh;
