# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

#SUBDIRS = 
#EXTRA_DIST =

.cpp.o:
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< 
	$(CXXCOMPILE) -c $< 


lib_LIBRARIES = libSiconosSharedLibrary.a
libSiconosSharedLibrary_a_SOURCES = SiconosSharedLibraryException.cpp SiconosSharedLibrary.cpp SiconosSharedLibraryException.h SiconosSharedLibrary.h 
include_HEADERS = SiconosSharedLibrary.h SiconosSharedLibraryException.h
