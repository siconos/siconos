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
	pfc_2D_solver.c \
	pfc_2D_cpg.c pfc_2D_latin.c pfc_2D_nlgs.c \
	lcp2dfc_2D.c dfc_2D2lcp.c \
	dfc_2D_solver.c \
	dfc_2D_latin.c \
	pr_solver.c \
	pr_nlgs.c pr_latin.c \
	dr_solver.c \
        dr_latin.c dr_nlgs.c \
	ql0001.f


libsolverpack_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

include_HEADERS = solverpack.h SiconosNumerics_lab.h # solverpack.h header to be added !!!
