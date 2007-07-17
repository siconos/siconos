# Siconos-Front-End version 2.1.1, Copyright INRIA 2005-2007.
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

include ../config/Make-rules

#include $(SCIDIR)/Makefile.incl



all     : Makelib
	make -f Makelib

Makelib : loader.sce

loader.sce :builder.sce 
	@echo running builder: BE PATIENT ...
	@$(SCIDIR)/util/scidem $(SCIDIR) builder.sce builder.dia

tests   : all
	make -f Makelib tests

clean   ::
	@echo Cleaning
	@$(RM) -f -r *.o *.lo *.la *.so *.dll *.sl *.a *.dia loader.sce libSiconosApiSci.c
	@$(RM) -f -r .libs f3cloader.sce f3cmake f3floader.sce f3fmake

distclean:: clean
	@$(RM) -f *.dia  *_gateway.c
	@$(RM) -f Makelib

message:
	@echo "------------------------------------------";
	@echo "At Scilab prompt, enter:";
	@echo "-->exec loader.sce";
	@echo "------------------------------------------";

