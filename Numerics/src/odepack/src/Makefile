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

OBJ = opkdmain.o  opkda1.o opkda2.o

all : build

build : $(OBJ)

lib : 
	$(ARCH) $(ARCHFLAGS) $(ODEPACKLIB) $(OBJ)
	$(RANLIB) $(ODEPACKLIB)
	$(CP) *.a $(NUMERICS_LIB)

clean:
	rm -f *.o *.a


%.o: %.f
	$(FORTRAN) $(FOPTS) -c -o $@ $(FFLAGS) $<
