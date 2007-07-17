# Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

include $(top_srcdir)/config/Make-rules

SUBDIRS = test
#EXTRA_DIST =

all : DefaultPlugin.so

DefaultPlugin.o: DefaultPlugin.cpp
	$(CXX) $(CXXFLAGS) -c DefaultPlugin.cpp

DefaultPlugin.so: DefaultPlugin.o
	$(CXX) $(CXXFLAGS) $(CXX_DYNLIB) $< -o $@

clean:
	$(RM) -f DefaultPlugin.o DefaultPlugin.so

pkgdata_DATA = DefaultPlugin.so


