# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUFFIXES = .cpp
.cpp.o:
	$(CXXCOMPILE) -c $<

bin_PROGRAMS = BALL2D
BALL2D_SOURCES = Ball2D.cpp 

BALL2D_LDADD = $(SICONOS_LDADD)
BALL2D_LDFLAGS = $(SICONOS_LDFLAGS)
