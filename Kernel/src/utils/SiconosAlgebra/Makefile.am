# Makefile.am --- SICONOS - SiconosAlgebra

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $< 

SUBDIRS = test 

noinst_LTLIBRARIES = libSiconosAlgebra.la
libSiconosAlgebra_la_SOURCES = SiconosVector.cpp SimpleVector.cpp CompositeVector.cpp SiconosMatrix.cpp \
				SiconosVector.h SimpleVector.h CompositeVector.h SiconosMatrix.h
include_HEADERS = SiconosVector.h SimpleVector.h CompositeVector.h SiconosMatrix.h 
libSiconosAlgebra_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
