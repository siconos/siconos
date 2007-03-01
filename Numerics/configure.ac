# Siconos-Numerics version 2.0.1, Copyright INRIA 2005-2006.
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
dnl configure script for SICONOS
dnl V. ACARY, 2005

AC_PREREQ(2.59)
AC_INIT([Siconos-Numerics],[2.0.1],[siconos.gforge.inria.fr])

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_TARGET([])
AC_CONFIG_SRCDIR([src/SiconosNumerics.h])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h]) 

# SO version of all siconos shared libraries
SICONOSNUMERICS_SO_CURRENT=1
SICONOSNUMERICS_SO_REVISION=2
SICONOSNUMERICS_SO_AGE=0
AC_SUBST(SICONOSNUMERICS_SO_CURRENT)
AC_SUBST(SICONOSNUMERICS_SO_REVISION)
AC_SUBST(SICONOSNUMERICS_SO_AGE)

# Checks for compilers.


m4_define([AC_PROG_INSURE],[
# Insure++ is a proriatary code checker/prober
AC_ARG_ENABLE(insure,
 [AC_HELP_STRING([--enable-insure],[Compile with Insure++])],
 [
  AC_CHECK_PROG(insure,insure,insure,no)
  AC_MSG_CHECKING([------------------------------------------------])
  if test "x$insure" = "xinsure"; then
   if test -z $CXX; then
    CXX=insure
   fi;
   if test -z $CC; then
    CC=insure
   fi;
  fi;
 ]
)
])  # AC_PROG_INSURE

dnl AC_PROG_INSURE

AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77

CFLAGS=""
CXXFLAGS=""
AC_MSG_CHECKING([g++ version])
CXX_VER=`$CXX -dumpversion`
#  CXX_VERSION=`$CXX -v 2>&1 | grep 'gcc version'`
case $CXX_VER in
     ( 1.*.* | 2.*.* )
     AC_MSG_ERROR(gcc 3 or superior is required.)
     exit 1;
     ;;
esac
AC_MSG_RESULT([ok ($CXX_VER)])   


# checks for libtool
AC_PROG_LIBTOOL


# Checks for programs.




# AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PATH_PROG(CP, cp, /bin/cp, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(MV, mv, /bin/mv, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(RM, rm, /bin/rm, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(MKDIR, mkdir, /bin/mkdir, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(ECHO, echo, /bin/echo, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(CAT, cat, /bin/cat, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(CHMOD, chmod, /bin/chmod, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(SED, sed, /bin/sed, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(AR, ar, /usr/bin/ar, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(DIFF, diff, /usr/bin/diff, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(AWK, awk, /usr/bin/awk, $PATH:/bin:/usr/bin:/usr/ucb/bin)
AC_PATH_PROG(DOXYGEN, doxygen, /usr/bin/doxygen, $PATH:/bin:/usr/bin:/usr/ucb/bin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h])
case $CXX_VER in
     ( 3.*.* )
     AC_CHECK_HEADERS([g2c.h])
     ;;
esac

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor gettimeofday])


###################################
## all the --enable-* argument help
###################################

AC_ARG_ENABLE(cppunit,[  --enable-cppunit              use cppunit to pass unitary tests, default=no ])
AC_ARG_ENABLE(debug,  [  --enable-debug              compiling with debugging symbols (-g), default=no ])
AC_ARG_ENABLE(optimisation,  [  --enable-optimisation              compiling with optimiation symbols (O2), default=no ])
AC_ARG_ENABLE(warnings,[  --enable-warnings            enable compilation warnings, default=yes])
AC_ARG_ENABLE(error-on-warning,[  --enable-erroronwarning     treat all compile warnings as errors, default=no])

#################################
## Checks for external libraries
#################################

dnl Check for Lapack
ACX_LAPACK([
  dnl action-if-found: 
  # This defines LAPACK_LIBS, BLAS_LIBS, and FLIBS
], [
  dnl action-if-not-found: 
  AC_MSG_ERROR([Blas/Lapack was not found. 
*** This means Lapack and matrix support cannot be compiled. 
*** This makes this library unusable. Please get blas and lapack 
*** installed. If you do have these installed, use the options 
*** --with-blas=<libname> or --with-lapack=<libname> and/or set 
*** the env variable LDFLAGS to include the appropriate linker 
*** flags.])
])

# Cppunit : Warning, by default we skip ccpunit
AC_MSG_CHECKING([cppunit])
if test "$enable_cppunit" = yes; then
AC_MSG_RESULT(option --enable-cppunit selected : able to run test suite)       
AM_PATH_CPPUNIT(1.10.2,,AC_MSG_ERROR([Cannot find  cppunit in a good version.]))
else
 AC_MSG_RESULT(option --enable-cppunit not selected : Warning : unable to run test suite)  
fi

AM_CONDITIONAL([CPPUNIT_COND], [test "$cppunit" = yes])


#################################
## COMPILE FLAGS CXXFLAGS CFLAGS FFLAGS
#################################

# Default Flag
if test "$enable_warnings" = no; then
CXXFLAGS=""
CFLAGS=""
FFLAGS=""
else
CXXFLAGS="-Wall"
CFLAGS="-Wall"
FFLAGS="-Wall"
fi

#Optional flag
if test "$enable_erroronwarning" = yes; then
   AC_MSG_RESULT(option --enable-erroronwarning selected : compiling with the -Werror option)
   CXXFLAGS="-Werror $CXXFLAGS"
   CFLAGS="-Werror $CFLAGS"
   FFLAGS="-Werror $FFLAGS"
fi


if test "$enable_debug" = yes; then
   AC_MSG_RESULT(option --enable-debug selected : compiling with the -g option)
   CXXFLAGS="-g -fbounds-check $CXXFLAGS"
   CFLAGS="-g $CFLAGS"
   FFLAGS="-g $FFLAGS"
fi

if test "$enable_optimisation" = yes; then
   AC_MSG_RESULT(option --enable-optimisation selected : compiling with the -03 option)
   CXXFLAGS="-O3 $CXXFLAGS"
   CFLAGS="-O3 $CFLAGS"
   FFLAGS="-O3 $FFLAGS"
fi


#Special Compiler flag for system type

AC_MSG_CHECKING([Special Compiler flag for system type])
case "$target" in
    *-mingw32*)
      ARCH="WIN32"
      ;;
    *-linux*)
      CFLAGS="$CFLAGS -fPIC"
      FFLAGS="$FFLAGS -fPIC"
      CXXFLAGS="$CXXFLAGS -fPIC"
      CXX_DYNLIB="-shared"
      ARCH="Linux"
      ;;
    *-apple-darwin*)
      # The Flag -fno-common is mandatory on MacOsX   
      FFLAGS="$FFLAGS -fno-common"  
      CXXFLAGS="$CXXFLAGS -flat_namespace"
      # This flag is used for LIBTOOL
      CXX_DYNLIB="-dynamiclib"
      #  if we want to have module and not dylib librairy
      #CXX_DYNLIB="-module"
      ARCH="Apple-Darwin"
      ;; 
    *)
      CXXFLAGS="$CXXFLAGS"
      CXX_DYNLIB="-shared"
      ARCH="OS_${target_os}"
      ;;
esac


CXXFLAGS="$CXXFLAGS "

AC_MSG_RESULT("Architecture $ARCH  Compilation flag CXXFLAGS= $CXXFLAGS CFLAGS= $CFLAGS FFLAGS= $FFLAGS  Linking flag CXX_DYNLIB = $CXX_DYNLIB")

dnl LDFLAGS="$CXX_DYNLIB"

AC_SUBST(CXX_DYNLIB)

#################################
## LINK FLAGS
#################################
#SICONOS_EXT_INCLUDES
SICONOS_EXT_INCLUDES=""
case "$target" in
    *-linux*)
      
      ;;
    *-apple-darwin*)
      SICONOS_EXT_INCLUDES="-I/sw/include/ $SICONOS_EXT_INCLUDES"
      ;;
esac
AC_MSG_RESULT("SICONOS_EXT_INCLUDES = $SICONOS_EXT_INCLUDES")
AC_SUBST(SICONOS_EXT_INCLUDES)

#SICONOS_EXT_LDFLAGS


SICONOS_EXT_LDFLAGS="$LAPACK_LIBS $BLAS_LIBS "

dnl case "$target" in
dnl     *-linux*)
dnl     AC_F77_LIBRARY_LDFLAGS	
dnl       ;;
dnl     *-apple-darwin*)
dnl     # the flag for linking c anf f77 has to be defined 
dnl     # the macros  AC_F77_LIBRARY_LDFLAGS yields an error of link between C and F77
dnl     FLIBS="" 
dnl     case $CXX_VER in
dnl      ( 3.*.* )
dnl      FLIBS="-L/sw/lib -lg2c"
dnl      ;;
dnl     esac	
dnl      ;;
dnl esac

SICONOS_EXT_LDFLAGS="$SICONOS_EXT_LDFLAGS $LIBS $FLIBS"
AC_SUBST(SICONOS_EXT_LDFLAGS)
AC_MSG_RESULT("SICONOS_EXT_LDFLAGS = $SICONOS_EXT_LDFLAGS")

#################################
## DOXYGEN 
#################################

DOXYGEN_CFG=$srcdir/config/DoxygenNumerics.cfg
AC_SUBST(DOXYGEN_CFG)


#################################
## OUTPUT FILES
#################################

# output files
AC_CONFIG_FILES([Makefile
		 src/Makefile
                 src/utils/Makefile
		 src/odepack/Makefile
		 src/odepack/test/Makefile
		 src/NSSpack/Makefile
		 src/NSSpack/test/Makefile])
AC_OUTPUT


echo "
 Configuration:
   Target:                              ${target}
   Architecture				${ARCH}
   C Compiler:                          ${CC}
   C compiler flags:                    ${CFLAGS}
   C++ Compiler:                        ${CXX}
   C++ compiler flags:                  ${CXXFLAGS}
   FORTRAN Compiler:                    ${F77}
   FORTRAN compiler flags:              ${FFLAGS}
   Preprocessor flags:                  ${CPPFLAGS}
   Linking flags			${LDFLAGS}
   Siconos includes:                    ${SICONOS_EXT_INCLUDES}
   Siconos ld flags:                    ${SICONOS_EXT_LDFLAGS}"
echo
echo Type \'make \' to build the SICONOS platform
echo Type \'make doc\' to generate the documentation.
echo Type \'make install\' to install the platfrom
echo Type \'make check\' to run the test
echo
echo Alternatively, you can \'make all\' to use the default configuration.
