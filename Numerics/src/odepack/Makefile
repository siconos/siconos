
#--------------------------------------------------
#                  SICONOS PROJECT
#--------------------------------------------------
 
#--------------------------------------------------
# Author(s)      : J-M BARBIER
# Creation Date  :
# Version        :
# Copyright      : LMGC 2004
#--------------------------------------------------
# Description    : Makefile for ODEPACK 
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
include $(NUMERICSPATH)/Makefile.inc
#  The machine (platform) identifier to append to the library names
#
#PLAT =  
#ODEPACKLIB    = odepack$(PLAT).a

OBJ = ./src/opkdmain.o  ./src/opkda1.o ./src/opkda2.o

help :
	@ echo "  "
	@ echo "+- NUMERICS/WP2 ------------------------------------------------------------+"
	@ echo "|  NUMERICS - LMGC 2004                                                   |" 
	@ echo "|                                                                          |"
	@ echo "|  usage :                                                                 |"
	@ echo "|       all             builds and link application                        |"
	@ echo "|       build           compiles and puts executable in bin directory      |"
	@ echo "|       inc             put the header files in the include directory      |"
	@ echo "|       lib             creates the library                                |"
	@ echo "|       clean           cleans *.o and executables files                   |"
	@ echo "|       buildtest       builds unitary tests                               |"
	@ echo "|       runtest         run unitary tests                                  |"
	@ echo "|       cleantest       cleans ojects and files generate by tests          |"
	@ echo "|                                                                          |"
	@ echo "+--------------------------------------------------------------------------+"
	@ echo "  "
	
all: inc build lib

build : 
	cd src; $(MAKE) build
	

lib : 
	cd src; $(MAKE) lib

inc :
	cd src;cp *.h $(NUMERICS_INC);cd ..

clean : 
	( cd src ; $(MAKE) clean )

 $(ODEPACKLIB): $(OBJ)
	$(ARCH) $(ARCHFLAGS) $(ODEPACKLIB) $(OBJ)
	$(RANLIB) $(ODEPACKLIB)

buildtest: 
	( cd test ; $(MAKE) buildtest )

cleantest:
	( cd test ; $(MAKE) cleantest )

runtest:
	( cd test ; $(MAKE) runtest )

#cppunit: 
#	( cd test ; $(MAKE) cppunit )

depend:

%.o : %.f
	$(FORTRAN) $(FOPTS) -c $< -o $@ 

