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

OBJETS_LCP=gcp_lcp.o gsnl_lcp.o latin_lcp.o solve_lcp.o minf.o maxf.o
#a venir lemke_lcp.o

#a venir OBJETS_QP=qp.o solve_qp.o
OBJETS_CFP= cfp_gsnl.o cfp_latin.o minf.o maxf.o solve_cfp.o cfp_gcp.o projc.o projf.o

OBJETS_CFD=gcp_lcp.o gsnl_lcp.o cfd_lcp.o cfd_latin.o lcp_cfd.o minf.o maxf.o solve_cfd.o 

OBJETS_RP=rp_latin.o rp_gsnl.o solve_rp.o maxf.o minf.o

OBJETS_RD=solve_rd.o maxf.o minf.o


all : build lib

build : $(OBJETS_LCP) $(OBJETS_RP) $(OBJETS_CFP) $(OBJETS_CFD) $(OBJETS_RD)

lib : NumericsSolvers_lcp.a NumericsSolvers_rp.a NumericsSolvers_cfp.a \
		NumericsSolvers_cfd.a NumericsSolvers_rd.a
	cp *.a $(NUMERICS_LIB)

NumericsSolvers_lcp.a: $(OBJETS_LCP)
	ar cr libNumericsSolvers_lcp.a $(OBJETS_LCP)

NumericsSolvers_rp.a: $(OBJETS_RP)
	ar cr libNumericsSolvers_rp.a $(OBJETS_RP)

NumericsSolvers_cfp.a: $(OBJETS_CFP)
	ar cr libNumericsSolvers_cfp.a $(OBJETS_CFP)

NumericsSolvers_cfd.a: $(OBJETS_CFD)
	ar cr libNumericsSolvers_cfd.a $(OBJETS_CFD)

NumericsSolvers_rd.a: $(OBJETS_RD)
	ar cr libNumericsSolvers_rd.a $(OBJETS_RD)
clean:
	rm -f *.o *.a

%.o: %.c
	$(CC) -c -o $@ $(CCFLAGS) $<
