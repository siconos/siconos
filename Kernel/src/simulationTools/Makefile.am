# Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
libSimulationTools_la_SOURCES =  UnitaryRelation.cpp Topology.cpp TimeDiscretisation.cpp Solver.cpp \
				Event.cpp TimeDiscretisationEvent.cpp NonSmoothEvent.cpp EventFactory.cpp EventsManager.cpp \
				OneStepNSProblem.cpp QP.cpp Relay.cpp FrictionContact.cpp FrictionContact2D.cpp FrictionContact3D.cpp LCP.cpp \
				OneStepIntegrator.cpp Lsodar.cpp Moreau.cpp Simulation.cpp TimeStepping.cpp EventDriven.cpp SparseBlockMatrix.cpp\
				UnitaryRelation.h UnitaryRelationsSet.h Topology.h TimeDiscretisation.h Solver.h \
				Event.h TimeDiscretisationEvent.h NonSmoothEvent.h EventFactory.h EventsManager.h \
				OneStepNSProblem.h QP.h Relay.h FrictionContact.h FrictionContact2D.h FrictionContact3D.h LCP.h \
				OneStepIntegrator.h Lsodar.h Moreau.h Simulation.h TimeStepping.h EventDriven.h \
				SimulationTools.h SimulationTypeDef.h SparseBlockMatrix.h
libSimulationTools_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
siconos_HEADERS = UnitaryRelation.h UnitaryRelationsSet.h Topology.h TimeDiscretisation.h Solver.h \
		  Event.h TimeDiscretisationEvent.h NonSmoothEvent.h EventFactory.h EventsManager.h \
		  OneStepNSProblem.h QP.h Relay.h FrictionContact.h FrictionContact2D.h FrictionContact3D.h LCP.h \
		  OneStepIntegrator.h Lsodar.h Moreau.h Simulation.h TimeStepping.h EventDriven.h \
		  SimulationTools.h SimulationTypeDef.h SparseBlockMatrix.h
