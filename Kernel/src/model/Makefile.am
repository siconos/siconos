# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUBDIRS = test

.cpp.o:
	$(CXXCOMPILE) -c $<

#lib_LIBRARIES = libModel.a
#libModel_a_SOURCES = Model.cpp KernelDefaultConfig.cpp Model.h KernelDefaultConfig.h

lib_LTLIBRARIES = libModel.la
libModel_la_SOURCES = Model.cpp KernelDefaultConfig.cpp Model.h KernelDefaultConfig.h
libModel_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)
include_HEADERS = Model.h KernelDefaultConfig.h SiconosConst.h
