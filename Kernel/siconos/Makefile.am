# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUFFIXES = .cpp
.cpp.o:
	$(CXXCOMPILE) -c $<

bin_PROGRAMS = SICONOS
SICONOS_SOURCES = main_siconos.cpp main_siconos.h
HEADER_LIST = main_siconos.h
