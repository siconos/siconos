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

# Checks for python.
AC_PYTHON_DEVEL
AC_MSG_CHECKING([Python Checking])
#ACX_CHECK_PYTHON([AC_SUBST(pythondir) AC_SUBST(PYTHON_DIR)], [AC_MSG_ERROR([Python problem])])
#echo $pythondir
#echo $PYTHON_SIR


# Checks for swig-python.
AC_PROG_SWIG(1.3.21)
SWIG_ENABLE_CXX
SWIG_PYTHON

# Checks for kernel & Numerics paths


#################################
## all the --with-* argument help
#################################
AC_ARG_WITH(localnumerics,	
	[  --with-localnumerics=path        indicates a local path of numerics lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])
AC_ARG_WITH(localkernel,	
	[  --with-localkernel=path        indicates a local path of kernel lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])

# Numerics
AC_MSG_CHECKING([SICONOS Numerics])
ACX_CHECK_NUMERICS([AC_SUBST(NUMERICS_INCLUDES) AC_SUBST(NUMERICS_LIBRARIES) AC_SUBST(NUMERICS_PATH)], [AC_MSG_ERROR([Cannot find  Siconos Numerics. Please install it or use the option --with-localnumerics.])])
echo $NUMERICS_PATH
# kernel
AC_MSG_CHECKING([SICONOS Kernel])
ACX_CHECK_KERNEL([AC_SUBST(KERNEL_INCLUDES) AC_SUBST(KERNEL_LIBRARIES) AC_SUBST(KERNEL_PATH)], [AC_MSG_ERROR([Cannot find  Siconos Kernel. Please install it or use the option --with-localkernel.])])
echo $KERNEL_PATH


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
		 config/Makefile
		 python/Makefile]
		 python/pySiconos/Makefile)
AC_OUTPUT
