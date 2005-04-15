#--------------------------------------------------
#                  SICONOS PROJECT
#--------------------------------------------------

#--------------------------------------------------
# Author(s)      : J-M BARBIER
# Creation Date  :
# Version        : 1.1
# Copyright      : LMGC 2004
#--------------------------------------------------
# Description    : general makefile for Siconos/Numerics
#--------------------------------------------------
# Options        : _ all             builds and link application
#                  _ build           builds the source files
#                  _ lib             creates the library and put them in the include directory
#                  _ inc             puts the header files in the include directory
#                  _ depend          generates the dependancy in Makefile.dep
#                  _ cleandep        cleans the dependencies
#                  _ clean           cleans *.o, *.a and executables files
#                  _ doc             generates Doxygen html documentation
#                  _ cleandoc        cleans the documentation
#                  _ buildtest       builds unitary tests
#                  _ runtest         runs unitary tests
#                  _ cleantest       cleans ojects and files generate by tests
#                  _ run             executes the main program
#                  _ cppunit         runs cppunit tests
#--------------------------------------------------
# Compiled files : none
#--------------------------------------------------
# Modifications  : _ Add command "depend"
#                  _ Add command "subdirs"
#                  _ correction on all the makefiles of Numerics for quality
#                  _
# Date           : 17/11/04
#--------------------------------------------------

#-----------------------------------------------
#NUMERICSPATH = $(PWD)
include $(NUMERICSPATH)/Makefile.inc

#--- PATHS ----------------------------------------
SUBDIRS = src
SUBDIRSLIBS = src

#--- TARGETS ---------------------------------------------

help :
	@ echo "  "
	@ echo "+- NUMERICS/WP2 ------------------------------------------------------------+"
	@ echo "|  NUMERICS - INRIA 2004                                                   |" 
	@ echo "|                                                                          |"
	@ echo "|  usage :                                                                 |"
	@ echo "|       all             builds and link application                        |"
	@ echo "|       build           compiles and puts executable in bin directory      |"
	@ echo "|       inc             put the header files in the include directory      |"
	@ echo "|       lib             creates the library                                |"
	@ echo "|       depend          generates the dependancy in Makefile.dep           |"
	@ echo "|       cleandep        cleans the dependencies                            |"
	@ echo "|       clean           cleans *.o, *.a and executables files              |"
	@ echo "|       doc             generates Doxygen html documentation               |"
	@ echo "|       cleandoc        cleans the documentation                           |"
	@ echo "|       buildtest       builds unitary tests                               |"
	@ echo "|       runtest         run unitary tests                                  |"
	@ echo "|       cleantest       cleans ojects and files generate by tests          |"
	@ echo "|       cppunit         compiles tests with cppunit                        |"
	@ echo "|       runcppunit      runs tests with cppunit                            |"
	@ echo "|       run             executes the main program                          |"
	@ echo "|                                                                          |"
	@ echo "+--------------------------------------------------------------------------+"
	@ echo "  "


.PHONY : $(SUBDIRS)

depend :
	$(RM) $(MAKEFILEDEP)
	@echo "####" $(MAKEFILEDEP) "####">> $(MAKEFILEDEP)
	@echo "#### Generated automatically by make depend on `date` ####">> $(MAKEFILEDEP)
	@echo "" >> $(MAKEFILEDEP)
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) depend); done

cleandep :
	$(RM) $(MAKEFILEDEP)
	touch $(MAKEFILEDEP)

clean : 
	$(RMDIR) $(NUMERICS_INC) $(NUMERICS_LIB) $(NUMERICS_BIN)
	$(RM) core* $(MAKEFILEDEP)
	touch $(MAKEFILEDEP)
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) clean); done

#------------------------------------------------------
all : inc build lib #subdirs lib

inc :
	$(MKDIR) $(NUMERICS_INC) #$(NUMERICS_LIB) #$(NUMERICS_BIN)
	for i in $(SUBDIRSLIBS); do (cd $$i; $(MAKE) inc); done

build : inc #ext
	(cd src; $(MAKE) build)
	@ echo "  "
	@ echo "%% build done %%"
	@ echo "  "

lib :	build
	$(MKDIR) $(NUMERICS_LIB)
	for i in $(SUBDIRSLIBS); do (cd $$i; $(MAKE) lib); done
	cd lib; $(MKDIR) tmp; cd tmp; $(ARCH) x ../$(ODEPACKLIB); $(ARCH) x ../$(SOLVERPACKLIB); \
	 gcc $(SHAREDFLAG) *.o -o $(LIBDYN_NUMERICS); cd ..; mv tmp/$(LIBDYN_NUMERICS) .; $(RMDIR) tmp
	@ echo "  "
	@ echo "%% lib done %%"
	@ echo "  "

#--------------------------------------------------
run :
	$(NUMERICS_EXE)

dist : clean cleantest
	(cd ..; tar cfvz save_NUMERICS.tar.gz Numerics/ )

#----- doc -------------------------------------------
doc : 
	$(MKDIR) doc
	cd doc; $(MKDIR) api
	$(DOXYGEN) $(DOXYGEN_CFG)


cleandoc :
	$(RMDIR) $(NUMERICS_DOC)


#--------- tests ----------------------------------
buildtest :
	@ echo  $(NUMERICSPATH)
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) buildtest); done
	@ echo "  "
	@ echo "     %% build test done %%"
	@ echo "  "

runtest :
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) -s runtest); done

cleantest :
	cd test; $(MAKE) cleantest
	for i in $(SUBDIRS); do (cd $$i; $(MAKE) cleantest); done

cppunit : buildtest
	cd test; $(MAKE) cppunit
	@ echo "  "
	@ echo "     %% build cppunit tests... done %%"
	@ echo "  "
	
runcppunit : cppunit
	cd test; $(MAKE) runcppunit
