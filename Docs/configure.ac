 # Siconos-Docs version 2.1.0, Copyright INRIA 2005-2006.
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
dnl configure script for SICONOS Docs

AC_PREREQ(2.59)
AC_INIT([Siconos-Docs],[2.1.0],[siconos.gforge.inria.fr])

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_TARGET([])
AC_CONFIG_SRCDIR([Makefile.am])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

#Checks libtool
AC_PROG_LIBTOOL

# Checks for programs.

# DOXYGEN
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([Siconos])

###################################
## all the --enable-* argument help
###################################

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


AC_CONFIG_FILES([Makefile
		 HomePage/config/Makefile 
		 User/GettingStarted/config/Makefile
                 User/InstallGuide/config/Makefile
		 User/UsersGuide/config/Makefile])

AC_OUTPUT
