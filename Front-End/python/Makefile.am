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
# Makefile.am --- SICONOS

SUBDIRS = pySiconos

include $(top_srcdir)/config/Make-rules
include $(KERNEL_PATH)/share/siconos-kernel/KernelConfig.sh
EXTRA_DIST = pySiconos.i 

all:swig_interface lib

swig_interface:pySiconos.i

pySiconos_wrap.cxx:pySiconos.i
	$(SWIG) $(SWIG_PYTHON_OPT) -I$(SICONOS_INCLUDES) pySiconos.i;
	sed "s/SWIG_append_msg/SWIG_append_errmsg/g" pySiconos_wrap.cxx > pySiconos_wrap.cxx.res ;
	mv pySiconos_wrap.cxx.res pySiconos_wrap.cxx;

pySiconos_wrap.o:pySiconos_wrap.cxx
	g++ -c $(CXXFLAGS) -I$(SICONOS_INCLUDES) $(LIBEXT_INCLUDES) pySiconos_wrap.cxx $(PYTHON_CPPFLAGS)

#wrap: libsiconos runswig pySiconos_wrap.o lib

#libsiconos:
#	./createLibSiconos.sh
#	g++ -g -shared -o libSiconos.so ./foo/*.o $(SICONOS_LDFLAGS) $(SICONOS_LDADD) $(LIBS)
#	rm -rf ./foo

#runswig: pySiconos.i
#	swig -python -c++ pySiconos.i;

#pySiconos_wrap.o: pySiconos_wrap.cxx
#	./correctif.sh;
#	g++ -c $(CXXFLAGS) $(INCLUDES) pySiconos_wrap.cxx -I/usr/include/python2.3;	

lib:_pySiconos.so

_pySiconos.so:pySiconos_wrap.o
	$(CXX) $(CXX_DYNLIB)  $(SICONOS_LDFLAGS)  pySiconos_wrap.o -o _pySiconos.so

#	g++ -shared $(SICONOS_LDFLAGS) $(SICONOS_LDADD) $(LIBS) libSiconos.so pySiconos_wrap.o -o _pySiconos.so 

clean:
	rm -f pySiconos_wrap.o _pySiconos.so pySiconos_wrap.cxx pySiconos.py pySiconos.pyc

install-data-am:
	test -z "$(PYTHON_DIR)" || $(mkdir_p) "$(PYTHON_DIR)"
	@INSTALL@ -m 0644 _pySiconos.so $(PYTHON_DIR)
	@INSTALL@ -m 0644 pySiconos.py $(PYTHON_DIR)
