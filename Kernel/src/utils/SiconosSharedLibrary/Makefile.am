# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

#SUBDIRS = 
#EXTRA_DIST =

.cpp.o:
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< 
	$(CXXCOMPILE) -c $< 


#lib_LIBRARIES = libSiconosSharedLibrary.a
#libSiconosSharedLibrary_a_SOURCES = SiconosSharedLibraryException.cpp SiconosSharedLibrary.cpp SiconosSharedLibraryException.h SiconosSharedLibrary.h 
lib_LTLIBRARIES = libSiconosSharedLibrary.la
libSiconosSharedLibrary_la_SOURCES = SiconosSharedLibraryException.cpp SiconosSharedLibrary.cpp SiconosSharedLibraryException.h SiconosSharedLibrary.h 
libSiconosSharedLibrary_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
include_HEADERS = SiconosSharedLibrary.h SiconosSharedLibraryException.h
