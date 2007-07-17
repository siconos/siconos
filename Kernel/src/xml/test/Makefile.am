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
# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

check_PROGRAMS = test1
TESTS = $(check_PROGRAMS)
test1_SOURCES = BoundaryConditionXMLTest.cpp LagrangianDSXMLTest.cpp OneStepIntegratorXMLTest.cpp SiconosMemoryXMLTest.cpp \
		DynamicalSystemXMLTest.cpp LinearDSXMLTest.cpp RelayNSLXMLTest.cpp InteractionXMLTest.cpp LinearTIRXMLTest.cpp \
		SiconosDOMTreeToolsTest.cpp TimeDiscretisationXMLTest.cpp TestMain.cpp
siconos_HEADERS = BoundaryConditionXMLTest.h LagrangianDSXMLTest.h OneStepIntegratorXMLTest.h SiconosMemoryXMLTest.h \
		DynamicalSystemXMLTest.h LinearDSXMLTest.h RelayNSLXMLTest.h InteractionXMLTest.h LinearTIRXMLTest.h \
		SiconosDOMTreeToolsTest.h TimeDiscretisationXMLTest.h
test1_INCLUDES = $(CPPUNIT_CFLAGS)
test1_LDFLAGS = $(CPPUNIT_LIBS) $(SICONOS_LIBS) $(SICONOS_EXT_LIBS)
