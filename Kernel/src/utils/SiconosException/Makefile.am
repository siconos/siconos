# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

#SUBDIRS =
#EXTRA_DIST =

SUFFIXES = .cpp
.cpp.o:
	$(CXXCOMPILE) -c $<

lib_LIBRARIES = libSiconosException.a
libSiconosException_a_SOURCES =  SiconosException.cpp RuntimeException.cpp XMLException.cpp 
include_HEADERS = RuntimeException.h    SiconosException.h    XMLException.h
libSiconosException_CXXFLAGS = -fPIC 
