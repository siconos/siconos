# Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
libSiconosXml_la_SOURCES = SiconosDOMTreeTools.cpp SiconosMemoryXML.cpp NonSmoothSolverXML.cpp DynamicalSystemXML.cpp FirstOrderNonLinearDSXML.cpp \
		FirstOrderLinearDSXML.cpp LagrangianDSXML.cpp LagrangianLinearTIDSXML.cpp NonSmoothLawXML.cpp NewtonImpactNSLXML.cpp \
		NewtonImpactFrictionNSLXML.cpp ComplementarityConditionNSLXML.cpp RelationXML.cpp FirstOrderRXML.cpp FirstOrderLinearRXML.cpp \
		LagrangianRXML.cpp LagrangianLinearRXML.cpp\
		InteractionXML.cpp NonSmoothDynamicalSystemXML.cpp SimulationXML.cpp TimeDiscretisationXML.cpp OneStepIntegratorXML.cpp \
		LsodarXML.cpp MoreauXML.cpp OneStepNSProblemXML.cpp QPXML.cpp RelayNSLXML.cpp\
		SiconosModelXML.cpp\
		XMLTagsName.h SiconosDOMTreeTools.h SiconosMemoryXML.h NonSmoothSolverXML.h DynamicalSystemXML.h FirstOrderNonLinearDSXML.h FirstOrderLinearDSXML.h \
		LagrangianDSXML.h LagrangianLinearTIDSXML.h NonSmoothLawXML.h NewtonImpactNSLXML.h NewtonImpactFrictionNSLXML.h \
		ComplementarityConditionNSLXML.h RelationXML.h FirstOrderRXML.h FirstOrderLinearRXML.h LagrangianRXML.h LagrangianLinearRXML.h\
		InteractionXML.h NonSmoothDynamicalSystemXML.h SimulationXML.h TimeDiscretisationXML.h OneStepIntegratorXML.h \
		LsodarXML.h MoreauXML.h OneStepNSProblemXML.h LCPXML.h FrictionContactXML.h QPXML.h RelayNSLXML.h \
		SiconosModelXML.h XmlTools.h

libSiconosXml_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
siconos_HEADERS = XMLTagsName.h SiconosDOMTreeTools.h SiconosMemoryXML.h NonSmoothSolverXML.h DynamicalSystemXML.h FirstOrderNonLinearDSXML.h FirstOrderLinearDSXML.h \
		LagrangianDSXML.h LagrangianLinearTIDSXML.h NonSmoothLawXML.h NewtonImpactNSLXML.h NewtonImpactFrictionNSLXML.h \
		ComplementarityConditionNSLXML.h RelationXML.h FirstOrderRXML.h FirstOrderLinearRXML.h LagrangianRXML.h LagrangianLinearRXML.h\
		InteractionXML.h NonSmoothDynamicalSystemXML.h SimulationXML.h TimeDiscretisationXML.h OneStepIntegratorXML.h \
		LsodarXML.h MoreauXML.h OneStepNSProblemXML.h LCPXML.h FrictionContactXML.h QPXML.h RelayNSLXML.h \
		SiconosModelXML.h XmlTools.h

