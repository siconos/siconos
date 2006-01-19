# Siconos version 1.0, Copyright INRIA 2005.  
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.	
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public license as published by
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
dnl JB CHARLETY, 2005

AC_PREREQ(2.57)
AC_INIT(SICONOS, 1.0, www.siconos.com)

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(SICONOS, 1.0)
AC_CONFIG_SRCDIR([src/model/Model.h])
AC_CONFIG_HEADER([config.h])



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
SICONOS_SO_CURRENT=0
SICONOS_SO_REVISION=1
SICONOS_SO_AGE=0
AC_SUBST(SICONOS_SO_CURRENT)
AC_SUBST(SICONOS_SO_REVISION)
AC_SUBST(SICONOS_SO_AGE)
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
]) # AC_PROG_INSURE

dnl  AC_PROG_INSURE

AC_PROG_CC
AC_PROG_CXX
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
AC_ARG_ENABLE(nana,[  --enable-nana              use nana to design by contract to pass unitary tests, default=no ])
AC_ARG_ENABLE(debug,  [  --enable-debug              compiling with debugging symbols (-g), default=no ])
AC_ARG_ENABLE(optimisation,  [  --enable-optimisation              compiling with optimiation symbols (O2), default=no ])
AC_ARG_ENABLE(warnings,[  --enable-warnings            enable compilation warnings, default=yes])
AC_ARG_ENABLE(error-on-warning,[  --enable-erroronwarning     treat all compile warnings as errors, default=no])

#################################
## all the --with-* argument help
#################################
AC_ARG_WITH(localnumerics,	
	[  --with-localnumerics=path        indicates a local path of numerics lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])
AC_ARG_WITH(localnana,			
	[  --with-localnana=path        indicates a local path of nana lib if the library is not installed in standard directories (/usr, /usr/local, ...)])
AC_ARG_WITH(localblas,			
	[  --with-localblas=path	indicates a local path of blas lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])
AC_ARG_WITH(locallapack,			  
	[  --with-locallapack=path      indicates a local path of lapack lib if the library is not installed in standard directories (/usr, /usr/local, ...)])
AC_ARG_WITH(locallapackpp,		
	[  --with-locallapackpp=path	indicates a local path of lapackpp lib if the library is not installed in standard directories (/usr, /usr/local, ...)])
AC_ARG_WITH(locallibxml,		
	[  --with-locallibxml=path	indicates a local path of libxml lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])
AC_ARG_WITH(localcppunit,		
	[  --with-localcppunit=path	indicates a local path of cppunit lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])

#################################
## Checks for external libraries
#################################

# required version of the libraries






# Blas
BLAS_VER="3.0."
AC_MSG_CHECKING([for Blas])
ACX_CHECK_BLAS([AC_SUBST(BLAS_LIBRARIES)], [AC_MSG_ERROR([Blas not found])])
# exit 1;

# Lapack
LAPACK_VER="3.0."
AC_MSG_CHECKING([for Lapack])
ACX_CHECK_LAPACK([AC_SUBST(LAPACK_LIBRARIES)], [AC_MSG_ERROR([Lapack not found])])
# exit 1;

# Numerics
AC_MSG_CHECKING([for Siconos/Numerics])
ACX_CHECK_NUMERICS([AC_SUBST(NUMERICS_INCLUDES) AC_SUBST(NUMERICS_LIBRARIES) AC_SUBST(NUMERICS_PATH)], [AC_MSG_ERROR([Cannot find  Numerics. Please install it or use the option --with-localnumerics.])])

# Lapack++
LAPACKPP_VER="1.8.0"	# lib....so number
LAPACKPP_REL_VER="2.4.2" # Release number
AC_MSG_CHECKING([for Lapack++])
ACX_CHECK_LAPACKPP([AC_SUBST(LAPACKPP_INCLUDES) AC_SUBST(LAPCKPP_LIBRARIES)], [AC_MSG_ERROR([Cannot find  Lapack++ in a good version. Please install Lapack++, version >= $LAPACKPP_VER  or use the option --with-locallapackpp. Warning: this version number corresponds to liblapackpp.so index. See INSTALL for corresponding release. ])])
CXX_LAPACKPP_FLAG=#"-DLA_COMPLEX_SUPPORT"  
AC_SUBST(CXX_LAPACKPP_FLAG)



# Libxml2
LIBXML_VER="2.6.16"
AC_MSG_CHECKING([Libxml2])
ACX_CHECK_LIBXML2([AC_SUBST(LIBXML2_INCLUDES) AC_SUBST(LIBXML2_LIBRARIES)], [AC_MSG_ERROR([Cannot find  Libxml2 in a good version. Please install Libxml2, version >= $LIBXML_VER  or use the option --with-locallibxml.])])

# Nana Warning, by default we skip ccpunit
AC_MSG_CHECKING([Nana])
if test "$enable_nana" = yes; then
ACX_CHECK_NANA([AC_SUBST(NANA_INCLUDES) AC_SUBST(NANA_LIBRARIES)], [AC_MSG_ERROR([Cannot find  Nana. Please install Nana  or use the option --with-localnana.])])
CXX_WITHOUTNANA_FLAG=""
else
 AC_MSG_RESULT(option --enable-nana not selected : Warning : do not use nana)  
CXX_WITHOUTNANA_FLAG="-DWITHOUT_NANA"
fi

# Cppunit : Warning, by default we skip ccpunit
CPPUNIT_VER="1.10.2"
AC_MSG_CHECKING([cppunit])
if test "$enable_cppunit" = yes; then
AC_MSG_RESULT(option --enable-cppunit selected : able to run test suite)       
ACX_CHECK_CPPUNIT([AC_SUBST(CPPUNIT_INCLUDES) AC_SUBST(CPPUNIT_LIBRARIES)], [AC_MSG_ERROR([Cannot find  cppunit in a good version. Please install cppunit, version >= $CPPUNIT_VER.])])
else
 AC_MSG_RESULT(option --enable-cppunit not selected : Warning : unable to run test suite)  
fi


#################################
## COMPILE FLAGS CXXFLAGS
#################################

# Default Flag
if test "$enable_warnings" = no; then
CXXFLAGS="-fPIC"
else
CXXFLAGS="-Wall -fPIC"
fi

#Optional flag
if test "$enable_erroronwarning" = yes; then
   AC_MSG_RESULT(option --enable-erroronwarning selected : compiling with the -Werror option)
   CXXFLAGS="-Werror $CXXFLAGS"
fi


if test "$enable_debug" = yes; then
   AC_MSG_RESULT(option --enable-debug selected : compiling with the -g option)
   CXXFLAGS="-g -fbounds-check $CXXFLAGS"
fi

if test "$enable_optimisation" = yes; then
   AC_MSG_RESULT(option --enable-optimisation selected : compiling with the -02 option)
   CXXFLAGS="-O2 $CXXFLAGS"
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


CXXFLAGS="$CXXFLAGS $CXX_WITHOUTNANA_FLAG $CXX_LAPACKPP_FLAG"

AC_MSG_RESULT("Architecture $ARCH  Compilation flag CXXFLAGS= $CXXFLAGS  Linking flag CXX_DYNLIB = $CXX_DYNLIB")


AC_SUBST(CXX_DYNLIB)

#################################
## LINK FLAGS
#################################

#SICONOS_EXT_INCLUDES
SICONOS_EXT_INCLUDES="$NUMERICS_INCLUDES $NANA_INCLUDES $LAPACKPP_INCLUDES $BLAS_INCLUDES $LAPACK_INCLUDES $LIBXML2_INCLUDES"
AC_SUBST(SICONOS_EXT_INCLUDES)
AC_MSG_RESULT("SICONOS_EXT_INCLUDES = $SICONOS_EXT_INCLUDES")
LIBG2C=""
#SICONOS_EXT_LDFLAGS
case $CXX_VER in
     ( 3.*.* )
     LIBG2C="-lg2c"
     ;;
esac  



#SICONOS_EXT_LDFLAGS
SICONOS_EXT_LDFLAGS="$NUMERICS_LIBRARIES $NANA_LIBRARIES $LAPACKPP_LIBRARIES $BLAS_LIBRARIES $LAPACK_LIBRARIES $LIBXML2_LIBRARIES"
case "$target" in
    *-linux*)
      SICONOS_EXT_LDFLAGS="$SICONOS_EXT_LDFLAGS $LIBG2C -ldl"
      ;;
    *-apple-darwin*)
      SICONOS_EXT_LDFLAGS="$SICONOS_EXT_LDFLAGS -L/sw/lib $LIBG2C "
      ;;
esac
AC_SUBST(SICONOS_EXT_LDFLAGS)
AC_MSG_RESULT("SICONOS_EXT_LDFLAGS = $SICONOS_EXT_LDFLAGS")

# LD_LIBRARY_PATH for siconos script
AC_SUBST(LOCAL_LIB_PATH) 

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
                 config/Makefile
		 config/xmlschema/Makefile
		 siconos/Makefile
                 src/Makefile
                 src/model/Makefile
                 src/model/test/Makefile
                 src/modelingTools/Makefile
		 src/modelingTools/test/Makefile
	         src/simulationTools/Makefile
                 src/plugin/Makefile
		 src/plugin_test/Makefile
                 src/utils/Makefile
		 src/utils/SiconosException/Makefile
                 src/utils/SiconosAlgebra/Makefile
                 src/utils/SiconosAlgebra/test/Makefile
                 src/utils/SiconosMemory/Makefile
                 src/utils/SiconosMemory/test/Makefile
                 src/utils/SiconosSharedLibrary/Makefile
                 src/utils/contract/Makefile
                 src/xml/Makefile
                 src/xml/test/Makefile])
AC_OUTPUT


echo "
 Configuration:
   Target:                              ${target}
   Architecture				${ARCH}
   C Compiler:                          ${CC}
   C compiler flags:                    ${CFLAGS}
   C++ Compiler:                        ${CXX}
   C++ compiler flags:                  ${CXXFLAGS}
   F77 Compiler:                        ${F77}
   F77 compiler flags:                  ${FFLAGS}
   Preprocessor flags:                  ${CPPFLAGS}
   Linking flags			${LDFLAGS}
   Includes:                            ${SICONOS_EXT_INCLUDES}
"
echo Type \'make \' to build the SICONOS platform
echo Run \'make doc\' to generate the documentation.
echo Type \'make install\' to install the platfrom
echo Type \'make check\' torun the test
echo
echo Alternatively, you can \'make all\' to use the default configuration.
