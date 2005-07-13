# Makefile.am --- SICONOS - Model

include $(top_srcdir)/config/Make-rules

SUBDIRS = test

.cpp.o:
	$(CXXCOMPILE) -c $<

#lib_LIBRARIES = libModel.a
#libModel_a_SOURCES = Model.cpp KernelDefaultConfig.cpp Model.h KernelDefaultConfig.h

noinst_LTLIBRARIES = libModel.la
libModel_la_SOURCES = Model.cpp KernelDefaultConfig.cpp Model.h KernelDefaultConfig.h
include_HEADERS = Model.h KernelDefaultConfig.h SiconosConst.h
