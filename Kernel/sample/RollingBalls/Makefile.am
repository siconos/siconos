# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUFFIXES = .cpp
.cpp.o:
	$(CXXCOMPILE) -c $<

bin_PROGRAMS = ROLLINGBALLS
ROLLINGBALLS_SOURCES = RollingBalls.cpp 

ROLLINGBALLS_LDADD = $(SICONOS_LDADD)
ROLLINGBALLS_LDFLAGS = $(SICONOS_LDFLAGS)
