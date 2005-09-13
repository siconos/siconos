# Makefile.am --- SICONOS - ModelingTools

include $(top_srcdir)/config/Make-rules

.c.o:
	$(CCCOMPILE) -c $<
.f.o:
	$(F77COMPILE) -c $<

SUBDIRS =

noinst_LTLIBRARIES = libsolverpack.la
libsolverpack_la_SOURCES = solverpack.h \
	projc.c projf.c minf.c maxf.c \
	lcp_solver.c lcp_solver_block.c \
	lcp_lexicolemke.c lcp_nlgs.c lcp_nsqp.c lcp_cpg.c \
	lcp_qp.c lcp_lemke.c lcp_latin.c lcp_newton_min.c\
	solve_cfp.c \
	cfp_gcp.c cfp_latin.c cfp_gsnl.c \
	lcp2cfd.c cfd2lcp.c \
	solve_rp.c \
	rp_nlgs.c rp_latin.c \
	solve_cfd.c \
	cfd_latin.c \
        rd_latin.c rd_nlgs.c solve_rd.c \
	ql0001.f


libsolverpack_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

include_HEADERS = solverpack.h SiconosNumerics_lab.h # solverpack.h header to be added !!!
