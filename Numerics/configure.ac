dnl configure script for SICONOS
dnl V. ACARY, 2005

AC_PREREQ(2.57)
AC_INIT(SICONOS, 1.0, www.siconos.com)

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(SICONOS/Numerics, 1.0)
AC_CONFIG_SRCDIR([src/solverpack/src/solve_lcp.c])
AC_CONFIG_HEADER([config.h])

# checks for libtool
AC_PROG_LIBTOOL

##############################################
# sh function to verify shared library version
##############################################

# function checkLibraryVersion ()
# $1 : absolute path of the library (ex. : /usr/lib/liblapack.so)
# $2 : minimum version (ex. : 3.0.3)

function checkLibraryVersion () { 
     lib=$1; 
     ver=$2;
     vers=`ls -l $lib | awk '{ print $NF}' |  awk 'BEGIN {FS = "."; } { printf "%d", ($3 * 1000 + $4) * 1000 + $5;}'`
     verReq=`expr $ver  | awk 'BEGIN {FS = "."; } { printf "%d", ($1 * 1000 + $2) * 1000 + $3;}'`
     if test -n "$vers" && test "$vers" -ge $verReq ; then
	return 1;
     else
	return 0;
     fi	     
}

# function checkLibraryVersionNumber ()
# $1 : version number ex. : 3.2.1)
# $2 : minimum version (ex. : 3.0.3)

function checkLibraryVersionNumber () { 
     vers=`expr $1 |  awk 'BEGIN {FS = "."; } { printf "%d", ($1 * 1000 + $2) * 1000 + $3;}'`
     verReq=`expr $2  | awk 'BEGIN {FS = "."; } { printf "%d", ($1 * 1000 + $2) * 1000 + $3;}'`
     if test -n "$vers" && test "$vers" -ge $verReq ; then
	return 1;
     else
	return 0;
     fi	     
}

# SO version of all siconos shared libraries
SICONOSNUMERICS_SO_CURRENT=0
SICONOSNUMERICS_SO_REVISION=1
SICONOSNUMERICS_SO_AGE=0
AC_SUBST(SICONOSNUMERICS_SO_CURRENT)
AC_SUBST(SICONOSNUMERICS_SO_REVISION)
AC_SUBST(SICONOSNUMERICS_SO_AGE)

# Checks for programs.
AC_PROG_CXX
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
AC_PROG_CC
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
## all the --with-* argument help
#################################
AC_ARG_WITH(localblas,			
	[  --with-localblas=path	indicates a local path of blas lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])
AC_ARG_WITH(locallapack,			  
	[  --with-locallapack=path      indicates a local path of lapack lib if the library is not installed in standard directories (/usr, /usr/local, ...)])


#################################
## Checks for external libraries
#################################

# required version of the libraries

LAPACK_VER="3.0.3"
BLAS_VER="3.0."
LIBXML_VER="2.6.17"
LAPACKPP_VER="1.6.0"
CPPUNIT_VER="1.10.2"

AC_MSG_CHECKING([Blas])
ACX_CHECK_BLAS([AC_SUBST(BLAS_LIBRARIES)], [AC_MSG_ERROR([Blas not found])])
# exit 1;

AC_MSG_CHECKING([Lapack])
ACX_CHECK_LAPACK([AC_SUBST(LAPACK_LIBRARIES)], [AC_MSG_ERROR([Lapack not found])])
# exit 1;


AC_MSG_CHECKING([cppunit])
$cppunit = "no"
if test "$enable_cppunit" = yes; then
AC_MSG_RESULT(option --enable-cppunit selected : able to run test suite)       
ACX_CHECK_CPPUNIT([AC_SUBST(CPPUNIT_INCLUDES) AC_SUBST(CPPUNIT_LIBRARIES) $cppunit = "yes"], [AC_MSG_ERROR([Cannot find  cppunit in a good version. Please install cppunit, version >= $CPPUNIT_VER.])])
else
 AC_MSG_RESULT(option --enable-cppunit not selected : Warning : unable to run test suite)  
fi


AM_CONDITIONAL([CPPUNIT_COND], [test "$cppunit" = yes])


#################################
## COMPILE FLAGS CXXFLAGS CFLAGS FFLAGS
#################################

# Default Flag
if test "$enable_warnings" = no; then
CXXFLAGS="-fPIC"
CFLAGS="-fPIC"
FFLAGS="-fPIC"
else
CXXFLAGS="-Wall -fPIC"
CFLAGS="-Wall -fPIC"
FFLAGS="-Wall -fPIC"
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
   CFLAGS="-g -fbounds-check $CFLAGS"
   FFLAGS="-g -fbounds-check $FFLAGS"
fi

if test "$enable_optimisation" = yes; then
   AC_MSG_RESULT(option --enable-optimisation selected : compiling with the -02 option)
   CXXFLAGS="-O2 $CXXFLAGS"
   CFLAGS="-O2 $CFLAGS"
   FFLAGS="-O2 $FFLAGS"
fi


#Special Compiler flag for system type

AC_MSG_CHECKING([Special Compiler flag for system type])
case "$target" in
    *-mingw32*)
      ARCH="WIN32"
      ;;
    *-linux*)
      CXXFLAGS="$CXXFLAGS"
      CXX_DYNLIB="-shared"
      ARCH="Linux"
      ;;
    *-apple-darwin*)
      CXXFLAGS="$CXXFLAGS -flat_namespace"
      CXX_DYNLIB="-dynamiclib -flat_namespace -single-module -undefined suppress"
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


AC_SUBST(CXX_DYNLIB)

#################################
## LINK FLAGS
#################################

#SICONOS_EXT_INCLUDES
SICONOS_EXT_INCLUDES="$BLAS_INCLUDES $LAPACK_INCLUDES"
AC_SUBST(SICONOS_EXT_INCLUDES)
AC_MSG_RESULT("SICONOS_EXT_INCLUDES = $SICONOS_EXT_INCLUDES")

#SICONOS_EXT_LDFLAGS
SICONOS_EXT_LDFLAGS="$BLAS_LIBRARIES $LAPACK_LIBRARIES"
case "$target" in
    *-linux*)
      SICONOS_EXT_LDFLAGS="-lg2c $SICONOS_EXT_LDFLAGS"
      ;;
    *-apple-darwin*)
      SICONOS_EXT_LDFLAGS="-L/sw/lib -lg2c $SICONOS_EXT_LDFLAGS"
      ;;
esac
AC_SUBST(SICONOS_EXT_LDFLAGS)
AC_MSG_RESULT("SICONOS_EXT_LDFLAGS = $SICONOS_EXT_LDFLAGS")






#################################
## DOXYGEN 
#################################

DOXYGEN_CFG=$srcdir/config/DoxygenSiconos.cfg
AC_SUBST(DOXYGEN_CFG)


#################################
## OUTPUT FILES
#################################

# output files
AC_CONFIG_FILES([Makefile
                 test/Makefile
                 src/Makefile
		 src/odepack/Makefile
		 src/odepack/test/Makefile
		 src/solverpack/Makefile
		 src/solverpack/src/Makefile
		 src/solverpack/test/Makefile])
AC_OUTPUT


