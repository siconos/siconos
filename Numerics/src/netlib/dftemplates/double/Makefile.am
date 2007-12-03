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

.f.o:
	$(F77COMPILE) -c $<

SUBDIRS =  

noinst_LTLIBRARIES = libdftemplates.la
libdftemplates_la_SOURCES =  Utils.f \
			MatVec.f PSolve.f \
			CG.f CGREVCOM.f \
			BiCG.f BiCGREVCOM.f \
			BiCGSTAB.f BiCGSTABREVCOM.f \
			CGS.f CGSREVCOM.f \
			Cheby.f ChebyREVCOM.f \
			GMRES.f GMRESREVCOM.f \
			Jacobi.f JacobiREVCOM.f \
			QMR.f QMRREVCOM.f \
			SOR.f SORREVCOM.f \
			STOPTEST2.f\
			GetBreak.f Matnorm.f

libdftemplates_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE)

numerics_HEADERS = dftemplates.h



check_PROGRAMS = test1 

TESTS = $(check_PROGRAMS)
include_HEADERS = 
test1_SOURCES = Tester.f dnsychk.f dspdchk.f dutils.f  Matgen.f
test1_INCLUDES = 
test1_LDFLAGS =  $(CPPUNIT_LIBRARIES)  $(DFTEMPLATES_LIB) $(SICONOS_EXT_LDFLAGS)
