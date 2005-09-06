# Makefile.am --- SICONOS - ModelingTools

include $(top_srcdir)/config/Make-rules

.c.o:
	$(CCCOMPILE) -c $<
.f.o:
	$(F77COMPILE) -c $<

SUBDIRS =

noinst_LTLIBRARIES = libsolverpack.la
libsolverpack_la_SOURCES = cfd_latin.c cfp_gcp.c cfp_latin.c gsnl_lcp.c lcp_cfd.c \
	lexicolemke_lcp.c minf.c  projf.c qpnonsym_lcp.c rp_gsnl.c solve_cfd.c solve_lcp.c\
	solve_rp.c cfd_lcp.c cfp_gsnl.c gcp_lcp.c latin_lcp.c lemke_lcp.c maxf.c\
        projc.c qp_lcp.c rd_latin.c rp_latin.c solve_cfp.c solve_rd.c \
	ql0001.f \
	SiconosNumerics.h SiconosNumerics_lab.h

libsolverpack_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

include_HEADERS = SiconosNumerics.h SiconosNumerics_lab.h # solverpack.h header to be added !!!
