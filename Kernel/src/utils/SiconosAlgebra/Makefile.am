# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $< 

SUBDIRS = test 

#lib_LIBRARIES = libSiconosAlgebra.a
#libSiconosAlgebra_a_SOURCES = SiconosVectorException.cpp SiconosMatrixException.cpp NewSiconosVector.cpp SimpleVector.cpp CompositeVector.cpp SiconosMatrix.cpp \
#				SiconosVectorException.h SiconosMatrixException.h NewSiconosVector.h SimpleVector.h CompositeVector.h SiconosMatrix.h
#include_HEADERS = SiconosVectorException.h SiconosMatrixException.h NewSiconosVector.h SimpleVector.h CompositeVector.h SiconosMatrix.h 


lib_LTLIBRARIES = libSiconosAlgebra.la
libSiconosAlgebra_la_SOURCES = SiconosVectorException.cpp SiconosMatrixException.cpp NewSiconosVector.cpp SimpleVector.cpp CompositeVector.cpp SiconosMatrix.cpp \
				SiconosVectorException.h SiconosMatrixException.h NewSiconosVector.h SimpleVector.h CompositeVector.h SiconosMatrix.h
include_HEADERS = SiconosVectorException.h SiconosMatrixException.h NewSiconosVector.h SimpleVector.h CompositeVector.h SiconosMatrix.h 
libSiconosAlgebra_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
