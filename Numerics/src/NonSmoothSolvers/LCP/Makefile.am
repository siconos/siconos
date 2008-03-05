# Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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

include $(top_srcdir)/config/Make-rules

# fmax
AM_CFLAGS=-std=c99

.c.o:
	$(CCCOMPILE) -c $<
.f.o:
	$(F77COMPILE) -c $<

##############################################################

noinst_LTLIBRARIES = libLCP.la

libLCP_la_SOURCES = lcp_cpg.c lcp_lexicolemke.c lcp_qp.c lcp_nsqp.c lcp_psor.c lcp_latin.c lcp_latin_w.c \
		    lcp_newton_FB.c lcp_path.c lcp_newton_min.c lcp_pgs.c lcp_rpgs.c \
		    filter_result_LCP.c lcp_compute_error.c \
		    LCP_driver.c lcp_solver_block.c \
		    LCP_Solvers.h LinearComplementarity_Problem.h

# lcp_solver_pred.c lcp_solver_block_pred_vec.c

libLCP_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

numerics_HEADERS = LCP_Solvers.h LinearComplementarity_Problem.h

