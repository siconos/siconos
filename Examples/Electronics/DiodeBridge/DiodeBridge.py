
# /* Siconos-sample , Copyright INRIA 2005-2011.
#  * Siconos is a program dedicated to modeling, simulation and control
#  * of non smooth dynamical systems.	
#  * Siconos is a free software; you can redistribute it and/or modify
#  * it under the terms of the GNU General Public License as published by
#  * the Free Software Foundation; either version 2 of the License, or
#  * (at your option) any later version.
#  * Siconos is distributed in the hope that it will be useful,
#  * but WITHOUT ANY WARRANTY; without even the implied warranty of
#  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  * GNU General Public License for more details.
#  *
#  * You should have received a copy of the GNU General Public License
#  * along with Siconos; if not, write to the Free Software
#  * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#  *
#  * Contact: Vincent ACARY vincent.acary@inrialpes.fr 
# */
# //-----------------------------------------------------------------------
# //
# //  DiodeBridge  : sample of an electrical circuit involving :
# //	- a linear dynamical system consisting of an LC oscillator (1 µF , 10 mH)
# //	- a non smooth system (a 1000 Ohm resistor supplied through a 4 diodes bridge) in parallel
# //	  with the oscillator
# //
# //  Expected behavior : 
# //	The initial state (Vc = 10 V , IL = 0) of the oscillator provides an initial energy.
# //	The period is 2 Pi sqrt(LC) ~ 0,628 ms.
# //      The non smooth system is a full wave rectifier :
# //	each phase (positive and negative) of the oscillation allows current to flow
# //	through the resistor in a constant direction, resulting in an energy loss :
# //	the oscillation damps.
# //
# //  State variables : 
# //	- the voltage across the capacitor (or inductor)
# //	- the current through the inductor
# //
# //  Since there is only one dynamical system, the interaction is defined by :
# //	- complementarity laws between diodes current and voltage. Depending on
# //        the diode position in the bridge, y stands for the reverse voltage across the diode 
# //	  or for the diode current (see figure in the template file) 
# //	- a linear time invariant relation between the state variables and
# //	  y and lambda (derived from Kirchhoff laws)
# //
# //-----------------------------------------------------------------------


t0 = 0.0
T = 5.0e-3       # Total simulation time
h_step = 1.0e-6  # Time step
Lvalue = 1e-2  # inductance
Cvalue = 1e-6   # capacitance
Rvalue = 1e3    # resistance 
Vinit = 10.0    # initial voltage
Modeltitle = "DiodeBridge"


from matplotlib.pyplot import subplot, title, plot, grid, show

from Siconos.Kernel import FirstOrderLinearDS
from numpy import array, eye, empty

init_state = array([1,0]) 

init_state[0] = Vinit

A = array([[0,-1./Cvalue],[1./Lvalue,0]])

LSDiodeBridge=FirstOrderLinearDS(init_state, A)
