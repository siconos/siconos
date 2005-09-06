# Makefile.am --- SICONOS - ModelingTools

include $(top_srcdir)/config/Make-rules

.cpp.o:
	$(CXXCOMPILE) -c $<

SUBDIRS = test 

noinst_LTLIBRARIES = libodepack.la
libodepack_la_SOURCES =  opkda1.f opkda2.f opkdmain.f
libodepack_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

include_HEADERS = # odepack.h header to be added !!!
