# Makefile.am --- SICONOS - SiconosException

include $(top_srcdir)/config/Make-rules

#SUBDIRS =
#EXTRA_DIST =

SUFFIXES = .cpp
.cpp.o:
	$(CXXCOMPILE) -c $<

noinst_LTLIBRARIES = libSiconosException.la
libSiconosException_la_SOURCES =  SiconosException.cpp RuntimeException.cpp XMLException.cpp SiconosMemoryException.cpp \
				  SiconosVectorException.cpp SiconosMatrixException.cpp SiconosSharedLibraryException.cpp \
				  SiconosException.h RuntimeException.h XMLException.h SiconosMemoryException.h \
			          SiconosVectorException.h SiconosMatrixException.h SiconosSharedLibraryException.h 
libSiconosException_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
include_HEADERS = SiconosException.h RuntimeException.h XMLException.h SiconosMemoryException.h SiconosVectorException.h SiconosMatrixException.h SiconosSharedLibraryException.h
