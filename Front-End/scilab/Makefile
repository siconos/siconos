# Siconos version 1.0, Copyright INRIA 2005.  
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


SICONOSPATH=/local_home/pissard/Workspace/siconos/distrib
NUMERICSPATH=/local_home/pissard/Workspace/siconos/distrib


include $(SICONOSPATH)/share/SICONOS/KernelConfig.sh

INC_PATH=$(SICONOSPATH)/include
CFLAGS=-g -fbounds-check -w -Wno-deprecated -fPIC -DLA_COMPLEX_SUPPORT -DWITHOUT_NANA -I$(INC_PATH) -I$(NUMERICSPATH)/include -I/usr/local//include/lapackpp/ `xml2-config --cflags`
DYNFLAGS=-shared

all: libsiconos.so sictest

libsiconos.so:siconos.cpp
	g++  -g  $(CFLAGS) -c siconos.cpp
	g++   -g $(DYNFLAGS) -o libsiconos.so siconos.o

sictest:sictest.c libsiconos.so
	gcc -g sictest.c -o sictest -L. -lsiconos $(SICONOS_LIBS) $(LIBEXT_LDFLAGS)

clean:
	rm libsiconos.so siconos.o sictest

