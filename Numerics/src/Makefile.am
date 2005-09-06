# Makefile.am --- SICONOS - src

SUBDIRS =  odepack solverpack utils
EXTRA_DIST = 

lib_LTLIBRARIES = libSiconosNumerics.la 
libSiconosNumerics_la_LDFLAGS= -version-info $(SICONOSNUMERICS_SO_CURRENT):$(SICONOSNUMERICS_SO_REVISION):$(SICONOSNUMERICS_SO_AGE) -no-undefined
libSiconosNumerics_la_SOURCES=
libSiconosNumerics_la_LIBADD = $(top_builddir)/src/odepack/libodepack.la $(top_builddir)/src/solverpack/src/libsolverpack.la
include_HEADERS = SiconosNumerics.h