# Makefile.am --- SICONOS

include $(top_srcdir)/config/Make-rules

SUFFIXES = .cpp
.cpp.o:
	$(CXXCOMPILE) -c $<

bin_PROGRAMS = DOUBLECONTACT
DOUBLECONTACT_SOURCES = DoubleContact.cpp 

DOUBLECONTACT_LDADD = $(SICONOS_LDADD)
DOUBLECONTACT_LDFLAGS = $(SICONOS_LDFLAGS)
