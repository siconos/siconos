# Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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

# fmax
AM_CFLAGS=-std=c99

.c.o:
	$(CCCOMPILE) -c $<
.f.o:
	$(F77COMPILE) -c $<

##############################################################

noinst_LTLIBRARIES = libMLCP.la

libMLCP_la_SOURCES = mlcp_path.c mlcp_psor.c mlcp_rpsor.c mlcp_pgs.c mlcp_rpgs.c \
		     mlcp_compute_error.c mlcp_filter_result.c \
		     mlcp_solver.c \
		     MLCP_Solvers.h

libMLCP_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

numerics_HEADERS = MLCP_Solvers.h
