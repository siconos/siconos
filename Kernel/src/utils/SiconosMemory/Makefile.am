# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUBDIRS = test

.cpp.o:
	$(CXXCOMPILE) -c $< 


lib_LIBRARIES = libSiconosMemory.a
libSiconosMemory_a_SOURCES = SiconosMemoryException.cpp SiconosMemory.cpp SiconosMemoryException.h SiconosMemory.h
include_HEADERS = SiconosMemoryException.h SiconosMemory.h
