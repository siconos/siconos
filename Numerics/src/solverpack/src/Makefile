#--------------------------------------------------
#                  SICONOS PROJECT
#--------------------------------------------------

#--------------------------------------------------
# Author(s)      : Jean-Michel Barbier
# Creation Date  :
# Version        : 1.1
# Copyright      : LMGC 2004
#--------------------------------------------------
# Description    : makefile for all the sources
#--------------------------------------------------
# Options        : _ all                      builds and link application
#                  _ build                    builds the source files
#                  _ lib                      creates the library and put them in the include directory
#                  _ clean                    cleans *.o and executables files
#                  _ NumericsSolvers_lcp.a    creates NumericsSolvers_lcp library
#                  _ NumericsSolvers_rp.a     creates NumericsSolvers_rp library
#                  _ NumericsSolvers_cfp.a    creates NumericsSolvers_cfp library
#                  _ NumericsSolvers_cfd.a    creates NumericsSolvers_cfd library
#--------------------------------------------------
# Compiled files :
#--------------------------------------------------
# Modifications  : _ correction on all the makefiles of Numerics for quality
#                  _ 
#                  _
# Date           : 17/11/2004
#--------------------------------------------------

#-----------------------------------------------
include $(NUMERICSPATH)/Makefile.inc

OBJETS_LCP= gcp_lcp.o gsnl_lcp.o latin_lcp.o lemke_lcp.o solve_lcp.o \
	ql0001.o  qp_lcp.o qpnonsym_lcp.o lexicolemke_lcp.o

OBJETS_CFP= cfp_gsnl.o cfp_latin.o solve_cfp.o cfp_gcp.o projc.o projf.o

OBJETS_CFD= cfd_lcp.o cfd_latin.o lcp_cfd.o gcp_lcp.o gsnl_lcp.o latin_lcp.o lemke_lcp.o solve_cfd.o

OBJETS_RP= rp_latin.o rp_gsnl.o solve_rp.o

OBJETS_RD= rd_latin.o solve_rd.o

OBJETS = $(OBJETS_LCP) $(OBJETS_RP) $(OBJETS_CFP) $(OBJETS_CFD) $(OBJETS_RD)

all : build lib

build : $(OBJETS_LCP) $(OBJETS_RP) $(OBJETS_CFP) $(OBJETS_CFD) $(OBJETS_RD)

lib : NumericsSolvers_lcp.a NumericsSolvers_rp.a NumericsSolvers_cfp.a \
		NumericsSolvers_cfd.a NumericsSolvers_rd.a libSolverpack.a
	cp *.a $(NUMERICS_LIB)
#	cp *.so $(NUMERICS_LIB)

NumericsSolvers_lcp.a: $(OBJETS_LCP)
	$(ARCH) cr libNumericsSolvers_lcp.a $(OBJETS_LCP)

NumericsSolvers_rp.a: $(OBJETS_RP)
	$(ARCH) cr libNumericsSolvers_rp.a $(OBJETS_RP)

NumericsSolvers_cfp.a: $(OBJETS_CFP)
	$(ARCH) cr libNumericsSolvers_cfp.a $(OBJETS_CFP)

NumericsSolvers_cfd.a: $(OBJETS_CFD)
	$(ARCH) cr libNumericsSolvers_cfd.a $(OBJETS_CFD)

NumericsSolvers_rd.a: $(OBJETS_RD)
	$(ARCH) cr libNumericsSolvers_rd.a $(OBJETS_RD)

libSolverpack.a:$(OBJETS)
	$(ARCH) cr libSolverpack.a $(OBJETS)
#	gcc $(SHAREDFLAG) $(OBJETS) -o libSolvers.so

clean:
	rm -f *.o *.a

%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $<
%.o: %.f
	$(FF) -c -o $@ $(FFLAGS) $<
