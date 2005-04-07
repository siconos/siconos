#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
dnl configure script for SICONOS front-end
dnl JB CHARLETY, 2005

AC_PREREQ(2.57)
AC_INIT(SICONOS_Front-end, 1.0, www.siconos.com)

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(SICONOS-Front-end, 1.0)
AC_CONFIG_SRCDIR([python/pySiconos.i])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB
AC_PYTHON_DEVEL
AC_PROG_SWIG(1.3.21)
SWIG_ENABLE_CXX
SWIG_PYTHON

# Checks for kernel & Numerics paths

# SICONOSPATH
AC_MSG_CHECKING([SICONOS Kernel Installation])
if test "$SICONOSPATH" != ""; then
	AC_MSG_RESULT([yes]);	
else	
	AC_MSG_ERROR(Cannot find SICONOS Kernel. Please install SICONOS Kernel or set SICONOSPATH environment variable.)
	exit 1;
fi

AC_MSG_CHECKING([SICONOS Numerics Installation])
if test "$NUMERICSPATH" != ""; then
	AC_MSG_RESULT([yes]);	
else	
	AC_MSG_ERROR(Cannot find SICONOS NUMERICS. Please install SICONOS Numerics or set NUMERICSPATH environment variable)
	exit 1;
fi

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
		 config/Makefile
		 python/Makefile])
AC_OUTPUT
