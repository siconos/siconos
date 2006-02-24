# Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
# Makefile.am --- SICONOS

SUBDIRS = xmlschema
EXTRA_DIST = 

confdir=$(pkgdatadir)
conf_DATA = KernelConfig.sh
KernelConfig.sh:KernelConfig.sh.in Makefile
	sed -e 's?\@SICONOS_EXT_INCLUDES\@?$(SICONOS_EXT_INCLUDES)?g' \
	    -e 's?\@SICONOS_EXT_LDFLAGS\@?$(SICONOS_EXT_LDFLAGS)?g' \
	    -e 's?\@CXX\@?$(CXX)?g' \
	    -e 's?\@CXXFLAGS\@?$(CXXFLAGS)?g' \
            -e 's?\@CXX_DYNLIB\@?$(CXX_DYNLIB)?g' \
	       < $(srcdir)/KernelConfig.sh.in > KernelConfig.sh.tmp \
	&& mv KernelConfig.sh.tmp KernelConfig.sh
