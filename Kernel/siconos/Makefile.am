# Makefile.am --- SICONOS

#include $(top_srcdir)/config/Make-rules

#SUFFIXES = .cpp
#.cpp.o:
#	$(CXXCOMPILE) -c $<

#bin_PROGRAMS = SICONOS
#SICONOS_SOURCES = main_siconos.cpp main_siconos.h
#HEADER_LIST = main_siconos.h

pkgdata_DATA = MakefilePluginSample
#MakefilePlugin MakefileSample

bin_SCRIPTS = test-script
test-script:test.in
	sed -e 's/whatever//' < $(srcdir)/test.in > test-script
	chmod +x test-script

