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
# Options        : _ all             builds and link application
#                  _ build           builds the source files
#                  _ lib             creates the library and put them in the include directory
#                  _ inc             puts the header files in the include directory
#                  _ clean           cleans *.o and executables files
#                  _ doc             generates Doxygen html documentation
#                  _ buildtest       builds unitary tests
#                  _ runtest         runs unitary tests
#                  _ cleantest       cleans ojects and files generate by tests
#--------------------------------------------------
# Compiled files :
#--------------------------------------------------
# Modifications  : _ correction on all the makefiles of Numerics for quality
#                  _ 
#
# Date           : 17/11/2004
#--------------------------------------------------

#-----------------------------------------------
include $(NUMERICSPATH)/Makefile.inc


help:
	@ echo "  "
	@ echo "+- NUMERICS/WP2 ------------------------------------------------------------+"
	@ echo "|  NUMERICS - LMGC 2004                                                    |" 
	@ echo "|                                                                          |"
	@ echo "|  usage :                                                                 |"
	@ echo "|       all             cleans and builds application                      |"
	@ echo "|       build           compiles source code                               |"
	@ echo "|       lib             create the libraries with the compiled source code |"
	@ echo "|       inc             puts the header files in the include directory     |"
	@ echo "|       doc             collect all the headers of the source file         |"
	@ echo "|       clean           cleans *.o and executables files                   |"
	@ echo "|       buildtest       builds unitary tests                               |"
	@ echo "|       runtest         run unitary tests                                  |"
	@ echo "|       cleantest       cleans ojects and files generate by tests          |"
#	@ echo "|       cppunit         run tests with cppunit                             |"
	@ echo "|                                                                          |"
	@ echo "+--------------------------------------------------------------------------+"
	@ echo "  "

depend:

#-------------------------------------
all: inc build lib

build :
	cd src; $(MAKE) build

lib : build
	cd src; $(MAKE) lib

inc :
	cd src;cp *.h $(NUMERICS_INC);cd ..

#-------------------------------------
doc:
	cd doc; rm -f headers.txt; head -80 ../src/*c > headers.txt

clean:
	cd src;	$(MAKE) clean

#-------------------------------------
buildtest : 
	cd test; $(MAKE) buildtest
	
runtest : 
	cd test; $(MAKE) runtest

cleantest:
	cd test ; $(MAKE) cleantest

#cppunit : 
#	cd test; $(MAKE) cppunit