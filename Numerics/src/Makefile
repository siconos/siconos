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
#                  _ buildtest       builds unitary tests
#                  _ runtest         runs unitary tests
#                  _ cleantest       cleans ojects and files generate by tests
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

#--- PATHS ----------------------------------------
SUBDIRS = odepack solverpack
SUBDIRSLIBS = odepack solverpack


#--------------------------------------------------

help :
	@ echo "  "
	@ echo "+- NUMERICS/WP2 ------------------------------------------------------------+"
	@ echo "|  NUMERICS - LMGC 2004                                                    |" 
	@ echo "|                                                                          |"
	@ echo "|  usage :                                                                 |"
	@ echo "|       all             cleans and builds application                      |"
	@ echo "|       build           compiles source code                               |"
	@ echo "|       lib             create the libraries with the compiled source code |"
	@ echo "|       inc             puts the header files in the include directory     |"
	@ echo "|       clean           cleans *.o and executables files                   |"
	@ echo "|       buildtest       builds unitary tests                               |"
	@ echo "|       runtest         run unitary tests                                  |"
	@ echo "|       cleantest       cleans ojects and files generate by tests          |"
	@ echo "|                                                                          |"
	@ echo "+--------------------------------------------------------------------------+"
	@ echo "  "

.PHONY: $(SUBDIRS)


depend : 
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) depend; cd ..;); done

#-------------------------------------

all : inc build lib

build : 
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) build); done

lib :
	for i in $(SUBDIRSLIBS); do (cd $$i; $(MAKE) lib); done

inc :
	for i in $(SUBDIRSLIBS); do (cd $$i; $(MAKE) inc); done

#-------------------------------------
clean :
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) clean); done

#-------------------------------------
buildtest :
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) buildtest); done

runtest :
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) runtest); done

#cppunit :
#	for i in $(SUBDIRS); do (cd $$i; $(MAKE) cppunit); done

cleantest :
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) cleantest); done
