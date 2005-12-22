# Siconos version 1.0, Copyright INRIA 2005.  
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
# Makefile.am --- SICONOS - simulationTools

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<


noinst_LTLIBRARIES = libSimulationTools.la
libSimulationTools_la_SOURCES = Solver.cpp LemkeSolver.cpp LexicoLemkeSolver.cpp QPSolver.cpp NSQPSolver.cpp NLGSSolver.cpp CPGSolver.cpp \
		LatinSolver.cpp Adams.cpp FrictionContact2D.cpp FrictionContact3D.cpp \
	  	EventDriven.cpp LCP.cpp Lsodar.cpp Moreau.cpp OneStepIntegrator.cpp OneStepNSProblem.cpp QP.cpp \
		Relay.cpp Strategy.cpp TimeDiscretisation.cpp TimeStepping.cpp \
		Solver.h LemkeSolver.h LexicoLemkeSolver.h QPSolver.h NSQPSolver.h NLGSSolver.h CPGSolver.h \
		LatinSolver.h Adams.h FrictionContact2D.h FrictionContact3D.h EventDriven.h LCP.h Lsodar.h Moreau.h OneStepIntegrator.h OneStepNSProblem.h QP.h \
		Relay.h Strategy.h TimeDiscretisation.h TimeStepping.h	
libSimulationTools_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
include_HEADERS = Solver.h LemkeSolver.h LexicoLemkeSolver.h QPSolver.h NSQPSolver.h NLGSSolver.h CPGSolver.h \
		LatinSolver.h Adams.h FrictionContact2D.h FrictionContact3D.h EventDriven.h LCP.h Lsodar.h Moreau.h OneStepIntegrator.h OneStepNSProblem.h QP.h \
		Relay.h Strategy.h TimeDiscretisation.h TimeStepping.h


