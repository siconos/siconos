# Makefile.am --- SICONOS - ModelingTools

include $(top_srcdir)/config/Make-rules

.c.o:
	$(CCCOMPILE) -c $<
.f.o:
	$(F77COMPILE) -c $<

SUBDIRS =

noinst_LTLIBRARIES = libsolverpack.la
libsolverpack_la_SOURCES = cfd_latin.c cfp_gcp.c cfp_latin.c lcp_cfd.c \
	lcp_lexicolemke.c lcp_nlgs.c lcp_nsqp.c lcp_solver.c lcp_cpg.c \
	lcp_qp.c lcp_lemke.c lcp_latin.c \
	minf.c  projf.c \
	solve_rp.c rp_gsnl.c solve_cfd.c \
	cfd_lcp.c cfp_gsnl.c maxf.c\
        projc.c rd_latin.c rp_latin.c solve_cfp.c solve_rd.c \
	ql0001.f \
	solverpack.h

libsolverpack_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

include_HEADERS = solverpack.h SiconosNumerics_lab.h # solverpack.h header to be added !!!
