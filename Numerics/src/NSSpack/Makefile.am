# Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2006.
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.	
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# Siconos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY vincent.acary@inrialpes.fr 
#	
# Makefile.am --- SICONOS - ModelingTools

include $(top_srcdir)/config/Make-rules

.c.o:
	$(CCCOMPILE) -c $<
.f.o:
	$(F77COMPILE) -c $<

SUBDIRS = test

### LCP #####################################################

LCP_SOURCES = lcp_solver.c lcp_solver_block_new.c \
	filter_result_LCP.c filter_result_LCP_block.c \
	lcp_compute_error.c $(LCP_SOLVERS)

LCP_SOLVERS = lcp_lexicolemke.c lcp_pgs.c lcp_nsqp.c lcp_cpg.c \
	lcp_qp.c lcp_latin.c lcp_latin_w.c lcp_newton_min.c \
	lcp_newton_FB.c lcp_sor.c

### PFC ######################################################

PFC_SOURCES = pfc_2D_solver.c pfc_3D_solver.c $(PFC_SOLVERS)

PFC_SOLVERS = pfc_2D_cpg.c pfc_2D_latin.c pfc_2D_nlgs.c \
	pfc_3D_cpg.c pfc_3D_nlgs.c \
	projc.c projf.c \
	pfc_2D_projc.c pfc_2D_projf.c \
	max_part.c min_part.c abs_part.c pos_part.c \
	pfc_3D_projc.c pfc_3D_projf.c

### DFC ######################################################

DFC_SOURCES = dfc_2D_solver.c $(DFC_SOLVERS)

DFC_SOLVERS = dfc_2D_latin.c   

### LCP <--> DFC 2D ##########################################

SWITCH = lcp2dfc_2D.c dfc_2D2lcp.c cond_2D2dfc_2D.c \
	dfc_2D2cond_2D.c 	


### PR #######################################################

PR_SOURCES = pr_solver.c $(PR_SOLVERS)

PR_SOLVERS = pr_nlgs.c pr_latin.c

### DR #######################################################

DR_SOURCES = dr_solver.c $(DR_SOLVERS)

DR_SOLVERS = dr_latin.c dr_nlgs.c

##############################################################

noinst_LTLIBRARIES = libNSSpack.la

libNSSpack_la_SOURCES = NSSpack.h \
	minf.c maxf.c \
	$(LCP_SOURCES) \
	$(PFC_SOURCES) \
	$(DFC_SOURCES) \
	$(SWITCH) \
	$(PR_SOURCES) \
	$(DR_SOURCES) \
	ql0001.f sortsn.f diffns.c  \
	freeSpBlMat.c


libNSSpack_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

numerics_HEADERS = NSSpack.h SiconosNumerics_lab.h
