#--------------------------------------------------
#                  SICONOS PROJECT
#--------------------------------------------------

#--------------------------------------------------
# Author(s)      :
# Creation Date  :
# Version        :
# Copyright      : INRIA 2004
#--------------------------------------------------

MAKEFILEDEP = $(SICONOSPATH)/Makefile.dep
include $(MAKEFILEDEP)

.SUFFIXES: .cpp .h .o .c

#--- C++ compiler and platform --------------------
#ARCH = Linux
GCC = g++
#GCC = insure

LAPACKPP_FLAGS = -DLA_COMPLEX_SUPPORT
CFLAGS = -g -fbounds-check -w -Wno-deprecated -fPIC $(LAPACKPP_FLAGS)

DYN_FLAGS = -shared

#--- LIBRARY -----------------------------------
LIB_CREATE = ar cr 
LIB_ADD = ar r 

#--- DIRECTORY ------------------------------------
SICONOS_DIR = $(SICONOSPATH)
SICONOS_LIB = $(SICONOS_DIR)/lib
SICONOS_INC = $(SICONOS_DIR)/include
SICONOS_BIN = $(SICONOS_DIR)/bin
SICONOS_EXT = $(SICONOS_DIR)/ext
SICONOS_DOC = $(SICONOS_DIR)/doc

#--- EXE ------------------------------------------
SICONOS_EXE = $(SICONOS_DIR)/siconos/SICONOS

#--- CLEANER --------------------------------------
RM = rm -f
RMDIR = rm -rf

#--- COPIER ---------------------------------------
CP = cp -r

#--- MOVE -----------------------------------------
MV = mv -f

#--- CREATE DIR -----------------------------------
MKDIR = mkdir -p

#--- MAKE -----------------------------------------
MAKE = make

#--- CPPUNIT --------------------------------------
#CPPUNIT_PATH = $(SICONOS_EXT)/cppunit
CPPUNIT_PATH = /home0/barbier/siconos/LIB_EXTERNALISATION/ext_EXT/cppunit-1.10.2/cppunit

CPPUNIT_INC = `$(CPPUNIT_PATH)/bin/cppunit-config --cflags` -w
CPPUNIT_LIB = `$(CPPUNIT_PATH)/bin/cppunit-config --libs`

#--- DOXYGEN --------------------------------------
DOXYGEN = doxygen
DOXYGEN_CFG = $(SICONOS_DIR)/config/DoxygenSiconos.cfg




#--- ODEPACK ------------------------------------
#ODEPACK_DIR = $(SICONOS_EXT)/ODEPACK
ODEPACK_DIR = /home0/barbier/siconos/LIB_EXTERNALISATION/ext_EXT/ODEPACK
ODEPACK_LIB = $(ODEPACK_DIR)/odepack.a

#--- NUMERICS ------------------------------------
NUMERICS = $(NUMERICSPATH)
NUMERICS_INC = -I$(NUMERICS)/include
NUMERICS_LIB = $(NUMERICS)/lib/*.a

#======================================
#--- NANA -----------------------------------------
NANA_DIR = /home0/barbier/siconos/LIB_EXTERNALISATION/ext_EXT/nana
NANA_INC = -I$(NANA_DIR)
NANA_LIB = $(NANA_DIR)/libnana.a
INCLUDENANA = $(NANA_INC)
WITHOUTNANA = -DWITHOUT_NANA
NANA_HOST =

#--- LIBF2C -----------------------------------
LIBF2C_DIR = /home0/barbier/siconos/LIB_EXTERNALISATION/ext_EXT/libf2c
LIBF2C_LIB = $(LIBF2C_DIR)/libf2c.a

#--- LIBXML -----------------------------------
#LIBXML_PATH = /home0/barbier/siconos/LIB_EXTERNALISATION/ext_EXT/libxml2
#LIBXML_INC = `$(LIBXML_PATH)/bin/xml2-config --cflags`
#LIBXML_LIB_PATH = -L $(LIBXML_PATH)/libxml2/lib
#LIBXML_LIB = -lxml2 -lz -lpthread -lm
#--- LIBXML -----------------------------------
#LIBXML_DIR = $(SICONOS_EXT)/libxml2
LIBXML_DIR = /home0/barbier/siconos/LIB_EXTERNALISATION/ext_EXT/libxml2
LIBXML_INC = `$(LIBXML_DIR)/bin/xml2-config --cflags`
LIBXML_LIB = `$(LIBXML_DIR)/bin/xml2-config --libs`

#--- LAPACK++ -----------------------------------
LAPACKPP_DIR = /home0/barbier/siconos/LIB_EXTERNALISATION/ext_EXT/lapackpp
LAPACKPP_INC = -I$(LAPACKPP_DIR)/include/lapackpp/
LAPACKPP_LIB = $(LAPACKPP_DIR)/lib/liblapack++.a $(LAPACKPP_DIR)/lib/libblas++.a $(LAPACKPP_DIR)/lib/liblamatrix++.a

#--- LAPACK ------------------------------------
LAPACK_DIR = /home0/barbier/siconos/LIB_EXTERNALISATION/ext_EXT/LAPACK
LAPACK_LIB = $(LAPACK_DIR)/lapack.a $(LAPACK_DIR)/blas.a

#======================================

#--- INCLUDE -----------------------------------
INCLUDE = -I$(SICONOS_INC) $(LIBXML_INC) $(LAPACKPP_INC) $(NANA_INC) $(NUMERICS_INC)
INCLUDE_EXT = -I$(SICONOS_INC) $(LIBXML_INC) $(LAPACKPP_INC) $(NANA_INC)
INCLUDE_DEP = -I$(SICONOS_DIR)/src/model -I$(SICONOS_DIR)/src/modelformalisation -I$(SICONOS_DIR)/src/modelstrategy \
			-I$(SICONOS_DIR)/src/utils -I$(SICONOS_DIR)/src/utils/SiconosSharedLibrary \
			-I$(SICONOS_DIR)/src/utils/SiconosMemory -I$(SICONOS_DIR)/src/utils/NewSiconosVector -I$(SICONOS_DIR)/src/xml \
			$(INCLUDE_EXT) $(CPPUNIT_INC) $(NUMERICS_INC)

#--- LIBS --------------------------------------
# original version 
LIBS = -ldl $(NANA_LIB) $(LAPACKPP_LIB) $(LAPACK_LIB) $(LIBF2C_LIB) \
		$(NUMERICS_LIB) $(LIBXML_LIB) # $(ODEPACK_LIB) 
#LIBS = -ldl $(NANA_LIB)   $(LAPACK_LIB) $(LAPACKPP_LIB)  $(ODEPACK_LIB)  $(LIBF2C_LIB)  $(LIBXML_LIB) # 

#---------------[ end of makefile.def ]-----------------------
