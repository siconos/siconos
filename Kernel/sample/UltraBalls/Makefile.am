# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUFFIXES = .cpp
.cpp.o:
	$(CXXCOMPILE) -c $<

bin_PROGRAMS = ULTRABALLS
ULTRABALLS_SOURCES = UltraBalls.cpp 

ULTRABALLS_LDADD = $(SICONOS_LDADD)
ULTRABALLS_LDFLAGS = $(SICONOS_LDFLAGS)
