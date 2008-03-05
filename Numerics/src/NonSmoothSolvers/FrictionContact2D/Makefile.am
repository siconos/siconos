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

noinst_LTLIBRARIES = libFrictionContact2D.la

libFrictionContact2D_la_SOURCES = cond_2D2dfc_2D.c dfc_2D2lcp.c lcp2dfc_2D.c pfc_2D_latin.c  pfc_2D_projc.c dfc_2D2cond_2D.c dfc_2D_latin.c\
				  dfc_2D_driver.c \
				  pfc_2D_cpg.c pfc_2D_nlgs.c pfc_2D_projf.c projc.c projf.c \
				  pfc_2D_driver.c \
				  FrictionContact2D_Solvers.h  
libFrictionContact2D_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

numerics_HEADERS = FrictionContact2D_Solvers.h 