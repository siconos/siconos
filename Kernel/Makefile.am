# Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

ACLOCAL AMFLAGS = -I macros

SUBDIRS = config src siconos global_tests man #src must be compiled before main program

# make dist with svn export
EXTRA_DIST = config src siconos global_tests test

doc:
	@if ! test -d $(top_srcdir)/../Docs; then \
	echo "Creation of ../Docs/htdocs/Doxygen."; \
	$(MKDIR) $(top_srcdir)/../Docs; \
	fi; \
	if ! test -d $(top_srcdir)/../Docs/htdocs; then \
	$(MKDIR) $(top_srcdir)/../Docs/htdocs; \
	fi; \
	if ! test -d $(top_srcdir)/../Docs/htdocs/Doxygen; then \
	$(MKDIR) $(top_srcdir)/../Docs/htdocs/Doxygen;  fi; \
	echo "Documentation files generated in ../Docs/htdocs/Doxygen ..."; 
	cd $(top_srcdir)/config; $(DOXYGEN) DoxygenSiconos.cfg

install-data-local: 
	  @if test -d $(top_srcdir)/doc; then \
	    if ! test -d @prefix@/doc; then \
	      $(MKDIR) @prefix@/doc; \
	    fi; \
	    if test -d @prefix@/doc/Siconos_Kernel; then \
	      $(RM) -rf @prefix@/doc/Siconos_Kernel; \
	    fi; \
	    $(CP) -r $(top_srcdir)/doc/ @prefix@/doc/Siconos_Kernel; \
	  fi; 
