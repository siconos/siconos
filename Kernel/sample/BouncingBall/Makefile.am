# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUFFIXES = .cpp
.cpp.o:
	$(CXXCOMPILE) -c $<

bin_PROGRAMS = BOUNCINGBALL
BOUNCINGBALL_SOURCES = BouncingBall.cpp 

BOUNCINGBALL_LDADD = $(SICONOS_LDADD)
BOUNCINGBALL_LDFLAGS = $(SICONOS_LDFLAGS)
