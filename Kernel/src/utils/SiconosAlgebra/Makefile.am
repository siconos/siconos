# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $< 


lib_LIBRARIES = libSiconosAlgebra.a
libSiconosAlgebra_a_SOURCES = SiconosVectorException.cpp SiconosMatrixException.cpp NewSiconosVector.cpp SimpleVector.cpp CompositeVector.cpp SiconosMatrix.cpp \
				SiconosVectorException.h SiconosMatrixException.h NewSiconosVector.h SimpleVector.h CompositeVector.h SiconosMatrix.h
include_HEADERS = SiconosVectorException.h SiconosMatrixException.h NewSiconosVector.h SimpleVector.h CompositeVector.h SiconosMatrix.h 

SUBDIRS = test 
