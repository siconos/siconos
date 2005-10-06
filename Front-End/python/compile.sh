#!/bin/sh

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

#
# Siconos lib
#
SRCP=../../../lib
LIBRARIES="$LIBP/lib_model.a $LIBP/lib_modelformalisation.a $LIBP/lib_modelstrategy.a $LIBP/lib_siconosMatrixVector.a $LIBP/lib_SiconosMemory.a $LIBP/lib_SiconosSharedLibrary.a $LIBP/lib_utils.a $LIBP/lib_xml.a"

TARGET=libSiconos

mkdir foo
cd foo
for lib in ${LIBRARIES}; do
    ar x $lib
done

rm ../${TARGET}.a
ar cr ../${TARGET}.a *
ranlib ../${TARGET}.a
g++ -shared -o ../${TARGET}.so *.o
cd ..
rm -rf foo

#
# Numeris lib
#
LIBP=../../../../Numerics/lib
LIBRARIES="$LIBP/NumericsSolvers_cfd.a $LIBP/NumericsSolvers_lcp.a $LIBP/NumericsSolvers_rp.a $LIBP/NumericsSolvers_cfp.a  $LIBP/NumericsSolvers_rd.a  $LIBP/odepack.a"
TARGET=libNumerics

rm ../${TARGET}.a
mkdir foo
cd foo
for lib in ${LIBRARIES}; do
    ar x $lib
done

ar cr ../${TARGET}.a *
ranlib ../${TARGET}.a
g++ -shared -o ../${TARGET}.so *.o
cd ..
rm -rf foo


g++ -g -fbounds-check -w -Wno-deprecated -fPIC -DLA_COMPLEX_SUPPORT -I/local_home/pissard/Workspace/SICONOS/include `xml2-config --cflags` -I/usr/local/include/lapackpp/ -I/usr/include -I/local_home/pissard/Workspace/Numerics/include -DWITHOUT_NANA -o .siconos main.o BouncingBall.o  libSiconos.so libNumerics.so  -ldl /usr/lib/libnana.a -llapack++ -lblas++ -llamatrix++ -llapack -lblas `xml2-config --libs` -lg2c -lm 
