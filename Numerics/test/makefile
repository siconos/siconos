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
#                  _ cppunit         launch cppunit tests
#                  _ cleantest       remove .o and executable files
#--------------------------------------------------
# Compiled files :
#--------------------------------------------------
# Modifications  : _ 
#                  _ 
#                  _
# Date           : 13/01/2005
#--------------------------------------------------

#-----------------------------------------------
include $(NUMERICSPATH)/Makefile.inc

#--- TESTEXE --------------------------------------
TESTEXE = ./TEST

#--- OBJ ------------------------------------------
TESTOBJ = TestNumerics.o
SRC = $(wildcard *.cpp)

all: cppunit runcppunit

cppunit : $(TESTEXE)

runcppunit : runtest


$(TESTEXE) : $(TESTOBJ)
	$(GCC) -g $(CPPUNIT_INC) $(INCLUDE) -o $(TESTEXE) TestMain.cpp $(TESTOBJ) $(CPPUNIT_LIB) $(LIBS)

$(TESTOBJ) : TestNumerics.cpp TestNumerics.h
	$(GCC) $(INCLUDE) $(CPPUNIT_INC) -c $(GCCFLAGS) TestNumerics.cpp

runtest :
	$(TESTEXE)
	
cleantest :
	rm -f $(TESTEXE) $(TESTOBJ) core.*
