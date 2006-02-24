# Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
# Makefile.am --- SICONOS - xml

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

#SUBDIRS = test


noinst_LTLIBRARIES = libSiconosXml.la
libSiconosXml_la_SOURCES = SolverXML.cpp BoundaryConditionXML.cpp FrictionContactXML.cpp \
		ComplementarityConditionNSLXML.cpp DSInputOutputXML.cpp DynamicalSystemXML.cpp \
		EqualityConstraintXML.cpp InteractionXML.cpp LagrangianDSIOXML.cpp LagrangianECXML.cpp LagrangianLinearDSIOXML.cpp \
		LagrangianLinearECXML.cpp LagrangianLinearRXML.cpp LagrangianDSXML.cpp LagrangianRXML.cpp \
		LagrangianLinearTIDSXML.cpp LCPXML.cpp LinearBCXML.cpp LinearDSIOXML.cpp \
		LinearDSXML.cpp LinearECXML.cpp LinearTIECXML.cpp LinearTIRXML.cpp\
		LsodarXML.cpp MoreauXML.cpp NewtonImpactLawNSLXML.cpp NewtonImpactFrictionNSLXML.cpp NLinearBCXML.cpp \
		NonSmoothLawXML.cpp NonSmoothDynamicalSystemXML.cpp OneStepIntegratorXML.cpp OneStepNSProblemXML.cpp \
		PeriodicBCXML.cpp QPXML.cpp RelationXML.cpp RelayNSLXML.cpp SiconosDOMTreeTools.cpp \
		SiconosMemoryXML.cpp SiconosModelXML.cpp StrategyXML.cpp TimeDiscretisationXML.cpp \
		SolverXML.h BoundaryConditionXML.h FrictionContactXML.h ComplementarityConditionNSLXML.h DSInputOutputXML.h \
		DynamicalSystemXML.h EqualityConstraintXML.h InteractionXML.h LagrangianDSIOXML.h LagrangianECXML.h LagrangianLinearDSIOXML.h \
		LagrangianLinearECXML.h LagrangianLinearRXML.h LagrangianDSXML.h LagrangianRXML.h \
		LagrangianLinearTIDSXML.h LCPXML.h LinearBCXML.h LinearDSIOXML.h LinearECXML.h \
		LinearDSXML.h LinearTIECXML.h LinearTIRXML.h LsodarXML.h MoreauXML.h \
		NewtonImpactLawNSLXML.h NewtonImpactFrictionNSLXML.h NLinearBCXML.h NonSmoothLawXML.h NonSmoothDynamicalSystemXML.h \
		OneStepIntegratorXML.h OneStepNSProblemXML.h PeriodicBCXML.h QPXML.h RelationXML.h \
		RelayNSLXML.h SiconosDOMTreeTools.h SiconosMemoryXML.h SiconosModelXML.h \
		StrategyXML.h TimeDiscretisationXML.h XMLTagsName.h

libSiconosXml_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
include_HEADERS = SolverXML.h BoundaryConditionXML.h FrictionContactXML.h ComplementarityConditionNSLXML.h \
		DSInputOutputXML.h DynamicalSystemXML.h EqualityConstraintXML.h InteractionXML.h \
		LagrangianDSIOXML.h LagrangianECXML.h LagrangianLinearDSIOXML.h LagrangianLinearECXML.h \
		LagrangianLinearRXML.h LagrangianDSXML.h \
		LagrangianRXML.h LagrangianLinearTIDSXML.h LCPXML.h LinearBCXML.h \
		LinearDSIOXML.h LinearECXML.h LinearDSXML.h LinearTIECXML.h LinearTIRXML.h \
		LsodarXML.h MoreauXML.h NewtonImpactLawNSLXML.h NewtonImpactFrictionNSLXML.h NLinearBCXML.h NonSmoothLawXML.h \
		NonSmoothDynamicalSystemXML.h OneStepIntegratorXML.h OneStepNSProblemXML.h PeriodicBCXML.h QPXML.h \
		RelationXML.h RelayNSLXML.h SiconosDOMTreeTools.h SiconosMemoryXML.h \
		SiconosModelXML.h StrategyXML.h TimeDiscretisationXML.h XMLTagsName.h
