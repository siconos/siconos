# Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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

SUBDIRS =  utils plugin modelingTools simulationTools controlTools model xml apiC
EXTRA_DIST = SiconosKernel.h

libSiconosKernel_la_LINK = $(LIBTOOL) --tag=CXX $(AM_LIBTOOLFLAGS) \
        $(LIBTOOLFLAGS) --mode=link $(CCLD) $(AM_CCFLAGS) \
        $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
lib_LTLIBRARIES= libSiconosKernel.la
libSiconosKernel_la_LDFLAGS= -version-info $(SICONOS_SO_CURRENT):$(SICONOS_SO_REVISION):$(SICONOS_SO_AGE) 
libSiconosKernel_la_SOURCES=
libSiconosKernel_la_LIBADD= $(top_builddir)/src/model/libModel.la \
	$(top_builddir)/src/simulationTools/libSimulationTools.la \
	$(top_builddir)/src/modelingTools/libModelingTools.la \
	$(top_builddir)/src/controlTools/libControlTools.la \
	$(top_builddir)/src/xml/libSiconosXml.la \
	$(top_builddir)/src/utils/SiconosMemory/libSiconosMemory.la \
	$(top_builddir)/src/utils/SiconosAlgebra/libSiconosAlgebra.la  \
	$(top_builddir)/src/utils/SiconosSharedLibrary/libSiconosSharedLibrary.la \
	$(top_builddir)/src/utils/SiconosException/libSiconosException.la\
	$(top_builddir)/src/apiC/libSiconosApiC.la

siconos_HEADERS = SiconosKernel.h


