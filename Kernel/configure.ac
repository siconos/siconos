# Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
# You should have received a copy of the GNU General public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY vincent.acary@inrialpes.fr 
#	
dnl configure script for SICONOS
dnl JB CHARLETY, 2005

AC_PREREQ(2.59)
AC_INIT([Siconos-Kernel],[2.1.0],[siconos.gforge.inria.fr])

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_TARGET([])
AC_CONFIG_SRCDIR([src/model/Model.h])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([macros])
# SO version of all siconos shared libraries
SICONOS_SO_CURRENT=2
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
   
dnl Uncomment the line above and add --enable-insure to compile with insure
dnl AC_PROG_INSURE

dnl Comment the following lines to compile with insure
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
     ( 4.1.*)
     CXXFLAGS="-ffriend-injection ";
     ;;   
esac
AC_MSG_RESULT([ok ($CXX_VER)])   

dnl End of comment for insure

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

dnl case "$target" in
dnl     *-apple-darwin*)
dnl      LDFLAGS="-L/sw/lib/"
dnl       ;; 
dnl esac
AC_CHECK_LIB(gmp, __gmpz_init, ,
       [AC_MSG_ERROR([GNU MP not found, see http://swox.com/gmp/. If the GNU MP is installed in a non standard directory DIR, set the variable LDFLAGS="-LDIR"])])

###################################
## all the --enable-* argument help
###################################

AC_ARG_ENABLE(cppunit,[  --enable-cppunit              use cppunit to pass unitary tests, default=no ])
AC_ARG_ENABLE(debug,  [  --enable-debug              compiling with debugging symbols (-g), default=no ])
AC_ARG_ENABLE(profiling,  [  --enable-profiling             compiling with profiling symbols (-pg), default=no ])
AC_ARG_ENABLE(optimisation,  [  --enable-optimisation              compiling with optimiation symbols (-O3), default=no ])
AC_ARG_ENABLE(warnings,[  --enable-warnings            enable compilation warnings, default=yes])
AC_ARG_ENABLE(error-on-warning,[  --enable-erroronwarning     treat all compile warnings as errors, default=no])

#################################
## all the --with-* argument help
#################################
AC_ARG_WITH(numerics,	
	[  --with-numerics=path        indicates a local path of numerics lib if the library is not installed in standard directories (/usr, /usr/local, ...) ])

dnl #################################
dnl ## Checks for external libraries
dnl #################################

# required version of the libraries

dnl dnl Check for Blas and Lapack
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

dnl Check for Atlas
AQ_CHECK_ATLAS
if test "$atlas_available" = "yes"; then
    AC_MSG_RESULT([Atlas OK])  
 else
    AC_MSG_ERROR([Cannot find  Atlas. Please install it or use the option  --with-atlas-libs=DIR  to add ATLAS library path.]) 
fi

case "$target" in
    *-linux*)
      ATLAS_LIBS="$atlas_libraries $atlas_lib"
      ATLAS_INCLUDES=""
      ;;
    *-apple-darwin*)
      #ATLAS_LIBS="-L/usr/local/lib/atlas/ -latlas -lcblas -lf77blas -llapack -latlas"
      ATLAS_LIBS="/usr/local/lib/atlas/libatlas.a /usr/local/lib/atlas/liblapack.a /usr/local/lib/atlas/libcblas.a /usr/local/lib/atlas/libf77blas.a"
      ATLAS_INCLUDES="-I/usr/local/include/atlas"
      ;; 
    *)
      ATLAS_LIBS="$atlas_libraries $atlas_lib"
      ATLAS_INCLUDES=""     
      ;;
esac


# Boost
AC_MSG_CHECKING([for Boost])
BOOST_VER="1.33.1"	# lib....so number
AX_BOOST([$BOOST_VER],,[AC_MSG_ERROR([Cannot find Boost in a good version. Please install Boost, version >= $BOOST_VER available at http://www.boost.org/.])])

# Libxml2
AC_MSG_CHECKING([for Libxml2])
AM_PATH_XML2(2.6.16,,AC_MSG_ERROR([
*** Libxml is required; Available at
http://www.libxml.org/.]))

# Cppunit : Warning, by default we skip ccpunit
AC_MSG_CHECKING([cppunit])
if test "$enable_cppunit" = yes; then
AC_MSG_RESULT(option --enable-cppunit selected : able to run test suite)       
AM_PATH_CPPUNIT(1.10,,AC_MSG_ERROR([Cannot find  cppunit in a good version.]))
else
 AC_MSG_RESULT(option --enable-cppunit not selected : Warning : unable to run test suite)  
fi

dnl Numerics
AC_MSG_CHECKING([for Siconos/Numerics])
ACX_CHECK_NUMERICS([AC_SUBST(NUMERICS_INCLUDES) AC_SUBST(NUMERICS_LIBRARIES) AC_SUBST(NUMERICS_PATH)], [AC_MSG_ERROR([Cannot find  Numerics. Please install it or use the option --with-numerics.])])


#################################
## COMPILE FLAGS CXXFLAGS
#################################
# Default Flag
if test "$enable_warnings" = no; then
CXXFLAGS="$CXXFLAGS"
else
CXXFLAGS="-Wall $CXXFLAGS"
fi

#Optional flag

if test "$enable_erroronwarning" = yes; then
   AC_MSG_RESULT(option --enable-erroronwarning selected : compiling with the -Werror option)
   CXXFLAGS="-Werror $CXXFLAGS"
fi

if test "$enable_profiling" = yes; then
	AC_MSG_RESULT(option --enable-profiling selected : compiling with the -pg option)
	CXXFLAGS="-pg -g -fbounds-check $CXXFLAGS"
else
	if  test "$enable_debug" = yes; then
	    AC_MSG_RESULT(option --enable-debug selected : compiling with the -g option)
	    CXXFLAGS="-g -fbounds-check $CXXFLAGS"
	else
		if test "$enable_optimisation" = yes; then
		   AC_MSG_RESULT(option --enable-optimisation selected : compiling with the -03 option)
		   CXXFLAGS="-DNDEBUG -O3 $CXXFLAGS"
		else
		CXXFLAGS="-DNDEBUG -O2 $CXXFLAGS"
		fi
	fi
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

CXXFLAGS="$CXXFLAGS"
#CXXFLAGS="$CXXFLAGS "

AC_MSG_RESULT("Architecture $ARCH  Compilation flag CXXFLAGS= $CXXFLAGS  Linking flag CXX_DYNLIB = $CXX_DYNLIB")


AC_SUBST(CXX_DYNLIB)

dnl #################################
dnl ## LINK FLAGS
dnl #################################

SICONOS_EXT_INCLUDES="$NUMERICS_INCLUDES $ATLAS_INCLUDES $XML_CPPFLAGS $BOOST_CPPFLAGS"
AC_SUBST(SICONOS_EXT_INCLUDES)


SICONOS_EXT_LIBS="$NUMERICS_LIBRARIES $FLIBS $ATLAS_LIBS $LAPACK_LIBS $XML_LIBS -ldl"
SICONOS_EXT_LIBS="$SICONOS_EXT_LIBS $LIBS " dnl  -L/usr/lib/atlas -llapack_atlas"

AC_SUBST(SICONOS_EXT_LIBS)


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
		 config/KernelConfig.sh
		 man/Makefile
		 man/siconos.1
		 siconos/MakefilePluginSample 
		 siconos/siconos
		 siconos/Makefile
                 src/Makefile
		 src/apiC/Makefile
                 src/model/Makefile
                 src/model/test/Makefile
                 src/modelingTools/Makefile
		 src/modelingTools/test/Makefile
	         src/simulationTools/Makefile
	         src/controlTools/Makefile
                 src/plugin/Makefile
		 src/plugin/test/Makefile
                 src/utils/Makefile
		 src/utils/SiconosException/Makefile
                 src/utils/SiconosAlgebra/Makefile
                 src/utils/SiconosAlgebra/test/Makefile
                 src/utils/SiconosMemory/Makefile
                 src/utils/SiconosMemory/test/Makefile
                 src/utils/SiconosSharedLibrary/Makefile
                 src/utils/contract/Makefile
		 src/utils/SiconosTools/Makefile
                 src/utils/boost/Makefile
		 src/xml/Makefile
                 src/xml/test/Makefile
		 global_tests/Makefile])
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
dnl   Linking flags			${LDFLAGS}
dnl   Linking flags(2)			${LIBS}
   Siconos external includes:           ${SICONOS_EXT_INCLUDES}
   Siconos lexternal libs :             ${SICONOS_EXT_LIBS}"
echo 
echo Type \'make \' to build the SICONOS platform
echo Type \'make doc\' to generate the documentation.
echo Type \'make install\' to install the platform
echo Type \'make check\' to run the test \(requires cppunit\)
echo 
echo Alternatively, you can \'make all\' to use the default configuration.

