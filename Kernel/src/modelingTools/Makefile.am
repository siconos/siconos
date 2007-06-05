# Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
libModelingTools_la_SOURCES = DynamicalSystem.cpp FirstOrderNonLinearDS.cpp FirstOrderLinearDS.cpp FirstOrderLinearTIDS.cpp LagrangianDS.cpp \
			      LagrangianLinearTIDS.cpp DynamicalSystemsSet.cpp Relation.cpp FirstOrderR.cpp FirstOrderLinearR.cpp FirstOrderLinearTIR.cpp LagrangianR.cpp \
			      LagrangianLinearR.cpp LagrangianScleronomousR.cpp LagrangianRheonomousR.cpp LagrangianCompliantR.cpp \
			      NonSmoothLaw.cpp NewtonImpactFrictionNSL.cpp NewtonImpactNSL.cpp RelayNSL.cpp ComplementarityConditionNSL.cpp \
			      Interaction.cpp InteractionsSet.cpp UnitaryRelation.cpp UnitaryRelationsSet.cpp Topology.cpp NonSmoothDynamicalSystem.cpp \
			      DynamicalSystem.h FirstOrderNonLinearDS.h FirstOrderLinearDS.h FirstOrderLinearTIDS.h LagrangianDS.h \
			      LagrangianLinearTIDS.h DynamicalSystemsSet.h Relation.h FirstOrderR.h FirstOrderLinearR.h FirstOrderLinearTIR.h LagrangianR.h LagrangianLinearR.h \
			      LagrangianScleronomousR.h LagrangianRheonomousR.h LagrangianCompliantR.h\
			      NonSmoothLaw.h NewtonImpactFrictionNSL.h NewtonImpactNSL.h RelayNSL.h ComplementarityConditionNSL.h \
			      Interaction.h InteractionsSet.h UnitaryRelation.h UnitaryRelationsSet.h Topology.h NonSmoothDynamicalSystem.h \
			      ModelingTools.h 
libModelingTools_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

siconos_HEADERS = DynamicalSystem.h FirstOrderNonLinearDS.h FirstOrderLinearDS.h FirstOrderLinearTIDS.h LagrangianDS.h \
		  LagrangianLinearTIDS.h DynamicalSystemsSet.h Relation.h FirstOrderR.h FirstOrderLinearR.h  FirstOrderLinearTIR.h LagrangianR.h LagrangianLinearR.h \
		  LagrangianScleronomousR.h LagrangianRheonomousR.h LagrangianCompliantR.h\
		  NonSmoothLaw.h NewtonImpactFrictionNSL.h NewtonImpactNSL.h RelayNSL.h ComplementarityConditionNSL.h \
		  Interaction.h InteractionsSet.h UnitaryRelation.h UnitaryRelationsSet.h Topology.h NonSmoothDynamicalSystem.h \
		  ModelingTools.h

