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
# Makefile.am --- SICONOS ModelingTools tests

include $(top_srcdir)/config/Make-rules
TESTS_ENVIRONMENT = LD_LIBRARY_PATH=$(top_srcdir)/src/plugin/:$(top_builddir)/src/plugin/test:${LD_LIBRARY_PATH}

check_PROGRAMS = test1
TESTS = $(check_PROGRAMS)
test1_SOURCES = FirstOrderRTest.cpp FirstOrderLinearRTest.cpp FirstOrderLinearTIRTest.cpp LagrangianScleronomousRTest.cpp LagrangianRheonomousRTest.cpp LagrangianCompliantRTest.cpp \
		LagrangianLinearRTest.cpp FirstOrderNonLinearDSTest.cpp FirstOrderLinearDSTest.cpp \
		LagrangianDSTest.cpp LagrangianLinearTIDSTest.cpp TestMain.cpp \
		FirstOrderRTest.h FirstOrderLinearRTest.h FirstOrderLinearTIRTest.h LagrangianScleronomousRTest.h  LagrangianRheonomousRTest.h LagrangianCompliantRTest.h \
	        LagrangianLinearRTest.h FirstOrderNonLinearDSTest.h FirstOrderLinearDSTest.h LagrangianDSTest.h LagrangianLinearTIDSTest.h
siconos_HEADERS =FirstOrderRTest.h FirstOrderLinearRTest.h FirstOrderLinearTIRTest.h LagrangianScleronomousRTest.h LagrangianRheonomousRTest.h LagrangianCompliantRTest.h \
		 LagrangianLinearRTest.h FirstOrderNonLinearDSTest.h FirstOrderLinearDSTest.h LagrangianDSTest.h LagrangianLinearTIDSTest.h
test1_INCLUDES = $(CPPUNIT_CFLAGS)
test1_LDFLAGS =  $(CPPUNIT_LIBS)  $(SICONOS_LIBS) $(SICONOS_EXT_LIBS)

