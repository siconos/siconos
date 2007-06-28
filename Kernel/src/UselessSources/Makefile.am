# Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
# Makefile.am --- SICONOS - ModelingTools

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

SUBDIRS = test 

noinst_LTLIBRARIES = libModelingTools.la
libModelingTools_la_SOURCES =  DynamicalSystemsSet.cpp InteractionsSet.cpp BoundaryCondition.cpp ComplementarityConditionNSL.cpp \
		DSInputOutput.cpp DynamicalSystem.cpp EqualityConstraint.cpp Interaction.cpp \
		Topology.cpp \
		LagrangianDSIO.cpp LagrangianEC.cpp LagrangianLinearDSIO.cpp LagrangianLinearEC.cpp \
		LagrangianLinearR.cpp LagrangianDS.cpp LagrangianR.cpp \
		LagrangianLinearTIDS.cpp LinearBC.cpp LinearDSIO.cpp LinearDS.cpp LinearTIDS.cpp\
		LinearEC.cpp LinearTIEC.cpp LinearTIR.cpp NewtonImpactNSL.cpp NewtonImpactFrictionNSL.cpp NLinearBC.cpp \
		NonSmoothLaw.cpp NonSmoothDynamicalSystem.cpp PeriodicBC.cpp Relation.cpp RelayNSL.cpp UnitaryRelation.cpp\
		UnitaryRelationsSet.cpp \
		DynamicalSystemsSet.h  InteractionsSet.h BoundaryCondition.h ComplementarityConditionNSL.h EqualityConstraint.h \
		DSInputOutput.h DynamicalSystem.h Interaction.h LagrangianDSIO.h LagrangianEC.h LagrangianLinearDSIO.h \
		LagrangianLinearEC.h LagrangianLinearR.h LagrangianDS.h LagrangianR.h LagrangianLinearTIDS.h \
		LinearBC.h LinearDSIO.h LinearEC.h LinearDS.h LinearTIDS.h LinearTIEC.h LinearTIR.h \
		NewtonImpactNSL.h NewtonImpactFrictionNSL.h NLinearBC.h NonSmoothLaw.h NonSmoothDynamicalSystem.h \
		PeriodicBC.h Relation.h RelayNSL.h UnitaryRelation.h UnitaryRelationsSet.h 
libModelingTools_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

siconos_HEADERS = DynamicalSystemsSet.h InteractionsSet.h BoundaryCondition.h ComplementarityConditionNSL.h EqualityConstraint.h \
		DSInputOutput.h DynamicalSystem.h Interaction.h \
		Topology.h \
		LagrangianDSIO.h LagrangianEC.h LagrangianLinearDSIO.h LagrangianLinearEC.h LagrangianLinearR.h \
		LagrangianDS.h LagrangianR.h LagrangianLinearTIDS.h LinearBC.h LinearDSIO.h \
		LinearEC.h LinearDS.h LinearTIDS.h LinearTIEC.h LinearTIR.h NewtonImpactNSL.h NewtonImpactFrictionNSL.h \
		NLinearBC.h NonSmoothLaw.h NonSmoothDynamicalSystem.h PeriodicBC.h Relation.h RelayNSL.h UnitaryRelation.h \
		UnitaryRelationsSet.h 
