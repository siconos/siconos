# Makefile.am --- SICONOS - SiconosMemory

include $(top_srcdir)/config/Make-rules

SUBDIRS = test

.cpp.o:
	$(CXXCOMPILE) -c $< 


noinst_LTLIBRARIES = libSiconosMemory.la
libSiconosMemory_la_SOURCES = SiconosMemoryException.cpp SiconosMemory.cpp SiconosMemoryException.h SiconosMemory.h
libSiconosMemory_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
include_HEADERS = SiconosMemoryException.h SiconosMemory.h
