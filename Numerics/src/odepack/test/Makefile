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
#                  _ testfunC        compiles test function : testFunC
#                  _ testfun-inC     compiles test function : testFun_inC
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

OBJ = DLSODE-test.o DLSODES-test.o DLSODA-test.o DLSODAR-test.o DLSODPK-test.o\
      DLSODKR-test.o DLSODI-test.o DLSOIBT-test.o DLSODIS-test.o 

EXE = DLSODE-test DLSODES-test DLSODA-test DLSODAR-test DLSODPK-test\
      DLSODKR-test  DLSODI-test DLSOIBT-test DLSODIS-test


all: build run

#cppunit :
	

buildtest : $(EXE) testfunC testfun-inC
	
runtest : $(EXE)
	for exefile in $(EXE); do \
	  ./$$exefile > $$exefile.log;	tail  $$exefile.log ; \
	done
	./test-funC
	./test-funC-inC

cleantest :
	( rm -f *.o *out *.log $(EXE) test-funC test-funC-inC)


$(EXE): $(OBJ)
	$(FORTRAN) $(FOPTS)   $@.o $(NUMERICS_LIB)/$(ODEPACKLIB) -o  $@

%.o : %.f
	$(FORTRAN) $(FOPTS) -c $< -o $@ 

%.o : %.c
	@ echo $(NUMERICSPATH)
	$(GCC) -c $(CFLAGS) $< $(LIBG2C_INC)

testfunC: funC.o DLSODE-test-funC.o
	$(FORTRAN) $(FOPTS)   DLSODE-test-funC.o  funC.o $(NUMERICS_LIB)/$(ODEPACKLIB) -o  test-funC -lstdc++

testfun-inC: funC.o test-funcC-inC.o
	$(GCC) $(CFLAGS) test-funcC-inC.o  funC.o $(NUMERICS_LIB)/$(ODEPACKLIB) -o  test-funC-inC $(LIBG2C_DIR) $(LIBG2C_LIB)
