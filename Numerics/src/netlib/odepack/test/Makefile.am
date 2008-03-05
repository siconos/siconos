# Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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

check_PROGRAMS = test1  test2 test3  test4  test5  test6  test7  test8  test9  test10 #test11

TESTS = $(check_PROGRAMS)
include_HEADERS = 
test1_SOURCES = DLSODES-test.f
test1_INCLUDES = $(CPPUNIT_INCLUDES)
test1_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test2_SOURCES = DLSODAR-test.f 
test2_INCLUDES = $(CPPUNIT_INCLUDES)
test2_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test3_SOURCES = DLSODI-test.f
test3_INCLUDES = $(CPPUNIT_INCLUDES)
test3_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test4_SOURCES = DLSODPK-test.f
test4_INCLUDES = $(CPPUNIT_INCLUDES)
test4_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test5_SOURCES = DLSODA-test.f
test5_INCLUDES = $(CPPUNIT_INCLUDES)
test5_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test6_SOURCES = DLSODE-test.f  
test6_INCLUDES = $(CPPUNIT_INCLUDES)
test6_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test7_SOURCES = DLSODIS-test.f   
test7_INCLUDES = $(CPPUNIT_INCLUDES)
test7_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test8_SOURCES =  DLSODKR-test.f  
test8_INCLUDES = $(CPPUNIT_INCLUDES)
test8_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test9_SOURCES = DLSOIBT-test.f 
test9_INCLUDES = $(CPPUNIT_INCLUDES)
test9_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

test10_SOURCES = test-funcC-inC.cpp funC.cpp
test10_INCLUDES = $(CPPUNIT_INCLUDES)
test10_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS)

# Change the linker fron g++ to g77 for this test.
#test11_SOURCES = DLSODE-test-funC.f funC.cpp
#test11_INCLUDES = $(CPPUNIT_INCLUDES)
#test11_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(ODEPACK_LIB) $(SICONOS_EXT_LDFLAGS) -lstdc++
