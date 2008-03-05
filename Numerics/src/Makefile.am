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
# Makefile.am --- SICONOS - src

include $(top_srcdir)/config/Make-rules

.c.o:
	$(CCCOMPILE) -c $<
.f.o:
	$(F77COMPILE) -c $<

SUBDIRS =  utils NonSmoothSolvers netlib
EXTRA_DIST = 

lib_LTLIBRARIES = libSiconosNumerics.la 
libSiconosNumerics_la_LDFLAGS=  $(CXX_DYNLIB) -version-info $(SICONOSNUMERICS_SO_CURRENT):$(SICONOSNUMERICS_SO_REVISION):$(SICONOSNUMERICS_SO_AGE) -release $(PACKAGE_VERSION)
libSiconosNumerics_la_SOURCES= dummy.cpp
libSiconosNumerics_la_LIBADD = 	$(top_builddir)/src/utils/libNumericsUtils.la \
				$(top_builddir)/src/NonSmoothSolvers/libNonSmoothSolvers.la \
				$(top_builddir)/src/netlib/odepack/libodepack.la \
			       	$(top_builddir)/src/netlib/dftemplates/double/libdftemplates.la 

numerics_HEADERS = SiconosNumerics.h

