# Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2006.
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

.c.o:
	$(CCCOMPILE) -c $<
.f.o:
	$(F77COMPILE) -c $<
.cpp.o:
	$(CXXCOMPILE) -c $<

include $(top_srcdir)/config/Make-rules

check_PROGRAMS = test1 test2 test3  test4  test5

#TESTS = $(check_PROGRAMS)
#XFAIL_TESTS = $(check_PROGRAMS)

TESTS = test1 test2 test3  test4  test5
XFAIL_TESTS =   test5


test1_SOURCES = main_lcp_new.c
include_HEADERS = 
test1_INCLUDES = $(CPPUNIT_INCLUDES)
test1_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(NSSPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test2_SOURCES =  main_dr.c
test2_INCLUDES = $(CPPUNIT_INCLUDES)
test2_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(NSSPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test3_SOURCES = main_pr.c 
test3_INCLUDES = $(CPPUNIT_INCLUDES)
test3_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(NSSPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test4_SOURCES = main_dfc_2D.c
test4_INCLUDES = $(CPPUNIT_INCLUDES)
test4_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(NSSPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test5_SOURCES = main_pfc_2D.c
test5_INCLUDES = $(CPPUNIT_INCLUDES)
test5_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(NSSPACK_LIB) $(SICONOS_EXT_LDFLAGS)

#test6_SOURCES = TestNumerics.cpp
#test6_INCLUDES = $(CPPUNIT_INCLUDES)
#test6_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(NSSPACK_LIB) $(SICONOS_EXT_LDFLAGS)
