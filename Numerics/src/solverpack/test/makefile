#--------------------------------------------------
#                  SICONOS PROJECT
#--------------------------------------------------

#--------------------------------------------------
# Author(s)      : Jean-Michel Barbier
# Creation Date  :
# Version        : 1.1
# Copyright      : LMGC 2004
#--------------------------------------------------
# Description    : makefile for all the sources
#--------------------------------------------------
# Options        : _ all             builds and runs the tests
#                  _ buildtest       builds the tests
#                  _ runtest         runs the tests
#                  _ cleantest       cleans built test files
#--------------------------------------------------
# Compiled files :
#--------------------------------------------------
# Modifications  : _ correction on all the makefiles of Numerics for quality
#                  _ 
#                  _
# Date           : 17/11/2004
#--------------------------------------------------

#-----------------------------------------------
include $(NUMERICSPATH)/Makefile.inc

#--- TESTEXE --------------------------------------
TESTEXE = ./TEST

#--- OBJ ------------------------------------------
OBJ =  $(NUMERICSPATH)/src/solverpack/src/*.o
TESTOBJ = main_rp.o main_cfd.o main_cfp.o main_lcp.o #TestNumerics.o
SRC = $(wildcard *.cpp)
LIB_LAPACK = $(NUMERICSPATH)/src/LAPACK/*.a

all: buildtest runtest

#cppunit : $(TESTEXE) TestMain.cpp runtest

#buildtest : $(TESTOBJ) #$(TESTEXE)
#	$(FF) -o $(TESTEXE) *.o $(NUMERICS_LIB)/*.a $(LIBS_IFC)

#$(TESTEXE) : $(TESTOBJ) #TestMain.cpp
#	$(GCC) $(INCLUDE) $(CPPUNIT_INC) -c $(GCCFLAGS) TestNumerics.cpp
#	$(GCC) -g $(CPPUNIT_INC) $(INCLUDE) -o $(TESTEXE) TestMain.cpp TestNumerics.o $(CPPUNIT_LIB) $(LIBS) $(NUMERICS_LIB)/*.a $(LIB_LAPACK)

buildtest :
#	$(GCC) $(INCLUDE) $(CPPUNIT_INC) -c $(GCCFLAGS) TestNumerics.cpp $(LIB_LAPACK)
	$(CC) $(INCLUDE) -c $(CCFLAGS) main_lcp.c
	$(FF) -o test_lcp main_lcp.o $(NUMERICS_LIB)/*.a $(LIBS_IFC)
	$(CC) $(INCLUDE) -c $(CCFLAGS) main_rp.c
	$(FF) -o test_rp main_rp.o $(NUMERICS_LIB)/*.a $(LIBS_IFC)
	$(CC) $(INCLUDE) -c $(CCFLAGS) main_cfd.c
	$(FF) -o test_cfd main_cfd.o $(NUMERICS_LIB)/*.a $(LIBS_IFC)
	$(CC) $(INCLUDE) -c $(CCFLAGS) main_cfp.c
	$(FF) -o test_cfp main_cfp.o $(NUMERICS_LIB)/*.a $(LIBS_IFC)

runtest :
	./test_cfp
	@ echo "## test_cfp done. ##"
	./test_cfd
	@ echo "## test_cfd done. ##"
	./test_rp
	@ echo "## test_rp done. ##"
	./test_lcp
	@ echo "## test_lcp done. ##"
	
cleantest :
	rm -f test_cfd test_cfp test_lcp test_rp $(TESTOBJ) core.*
