 # Siconos-Front-End version 3.0.0, Copyright INRIA 2005-2008.
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
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
dnl configure script for SICONOS front-end
dnl JB CHARLETY, 2005

AC_PREREQ(2.59)
AC_INIT([Siconos-Front-end],[3.0.0],[siconos.gforge.inria.fr])

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_TARGET([])
AC_CONFIG_SRCDIR([python/pySiconos.i])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

#Checks libtool
AC_PROG_LIBTOOL

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB

###################################
## all the --enable-* argument help
###################################

AC_ARG_ENABLE(scilab,[  --enable-scilab             compiling toolbox for scilab, default=no ])

# Checks for python.
AC_PYTHON_DEVEL
AC_MSG_CHECKING([Python Checking])
#ACX_CHECK_PYTHON([AC_SUBST(pythondir) AC_SUBST(PYTHON_DIR)], [AC_MSG_ERROR([Python problem])])
#echo $pythondir
#echo $PYTHON_SIR

dnl make /usr/local the default for the installation

AC_PREFIX_DEFAULT(/usr/local)
local_installation="no"
if test "x$prefix" = "xNONE"; then
  prefix=$ac_default_prefix
  ac_configure_args="$ac_configure_args --prefix $prefix"
 else
 local_installation="yes"
fi 
AC_SUBST(prefix)
echo "local_installation" $local_installation

# Checks for swig-python.
AC_PROG_SWIG(1.3.21)
SWIG_ENABLE_CXX
SWIG_PYTHON

# Scilab : Warning, by default we skip scilab
AC_MSG_CHECKING([scilab])
if test "$enable_scilab" = yes; then
# test if SCI is defined. If not,
# test if scilab is in path, then get SCI from scilab -nw.
if test -z $SCI; then
   AC_CHECK_PROG([has_scilab],[scilab],[yes],[no], , )
   if test x$has_scilab = xno; then
      AC_MSG_ERROR([[neither Scilab was found in your PATH nor a SCI variable has been defined]])
   else
	cmd='F=mopen("path.incl","w");
        mfprintf(F,SCI);
        mclose(F);'
      echo "$cmd" | scilab -nw >/dev/null
      SCIDIR=`cat path.incl` 
      echo "SCIDIR=" $SCIDIR
      SCILAB=scilab
      AC_SUBST(SCIDIR)
      AC_SUBST(SCILAB)
   fi
else
   SCIDIR=$SCI
   echo "SCIDIR=" $SCIDIR
   AC_SUBST(SCIDIR)
   AC_SUBST(SCILAB)
fi
else
 AC_MSG_RESULT(option --enable-scilab not selected : Warning : unable to compile scilab toolbox)  
fi

# Checks for kernel & Numerics paths



#################################
## all the --with-* argument help
#################################
AC_ARG_WITH(numerics,	
	[  --with-numerics=path        indicates a local path of numerics lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])
AC_ARG_WITH(kernel,	
	[  --with-kernel=path        indicates a local path of kernel lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])

# Numerics
AC_MSG_CHECKING([SICONOS Numerics])
ACX_CHECK_NUMERICS([AC_SUBST(NUMERICS_INCLUDES) AC_SUBST(NUMERICS_LIBRARIES) AC_SUBST(NUMERICS_PATH)], [AC_MSG_ERROR([Cannot find  Siconos Numerics. Please install it or set correctly the option --with-numerics.])])
echo "NUMERICS_PATH=" $NUMERICS_PATH
# kernel
AC_MSG_CHECKING([SICONOS Kernel])
ACX_CHECK_KERNEL([AC_SUBST(KERNEL_INCLUDES) AC_SUBST(KERNEL_LIBRARIES) AC_SUBST(KERNEL_PATH)], [AC_MSG_ERROR([Cannot find  Siconos Kernel. Please install it or set correctly the option --with-kernel.])])
echo "KERNEL_PATH=" $KERNEL_PATH

#PYTHONPATH


PYTHON_DIR=$PYTHON_SITE_PKG

if test "$local_installation" = "yes" ; then
   PYTHON_DIR=$prefix/site-package
   PYTHONPATH=$prefix/site-package
fi
echo "PYTHON_DIR=" $PYTHON_DIR
echo "PYTHONPATH=" $PYTHONPATH
AC_SUBST(PYTHON_DIR)
AC_SUBST(PYTHONPATH)

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
		 config/Makefile
		 python/Makefile
		 python/pySiconos/pySiconos
		 python/pySiconos/Makefile
		 python/pySiconos/man/Makefile
		 python/pySiconos/man/pySiconos.1
		 scilab/Makefile
		 scilab/loaderSiconos.sce
		 scilab/env.sh]
		 )
AC_OUTPUT
