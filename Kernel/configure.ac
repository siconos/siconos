dnl configure script for SICONOS
dnl JB CHARLETY, 2005


AC_PREREQ(2.57)
AC_INIT(SICONOS, 1.0, www.siconos.com)
AM_INIT_AUTOMAKE(SICONOS, 1.0)
AC_CONFIG_SRCDIR([src/model/Model.h])
AC_CONFIG_HEADER([config.h])

dnl Checking host/target/build systems, for make, install etc.
AC_CANONICAL_SYSTEM

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
AC_PROG_RANLIB
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

AC_ARG_ENABLE(cppunit,[  --enable-cppunit              use cppunit to pass unitary tests ])

#################################
## all the --with-* argument help
#################################

AC_ARG_WITH(localnana,			
	[  --with-localnana=path        indicates a local path of nana lib if the library is not installed in /usr or /usr/local ])
AC_ARG_WITH(localblas,			
	[  --with-localblas=path	indicates a local path of blas lib if the library is not installed in /usr or /usr/local])
AC_ARG_WITH(locallapack,			  
	[  --with-locallapack=path      indicates a local path of lapack lib if the library is not installed in /usr or /usr/local])
AC_ARG_WITH(locallapackpp,		
	[  --with-locallapackpp=path	indicates a local path of lapackpp lib if the library is not installed in /usr or /usr/local])
AC_ARG_WITH(locallibxml,		
	[  --with-locallibxml=path	indicates a local path of libxml lib if the library is not installed in /usr or /usr/local])

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

# Checks for external libraries

LAPACK_VER="3.0.3"
BLAS_VER="3.0."
LIBXML_VER="2.6.17"
LAPACKPP_VER="1.5.1"
CPPUNIT_VER="1.10.2"

#AC_MSG_CHECKING([Blas])
ACX_CHECK_BLAS([AC_SUBST(BLAS_LIBRARIES)], [AC_MSG_ERROR([Blas not found])])
# exit 1;

# Numerics
AC_MSG_CHECKING([SICONOS Numerics - headers])
for ac_dir in		\
    /usr/include	\
    /usr/local/include	\
   $NUMERICSPATH/include	\
    ; \
do  if test -r "$ac_dir/SiconosNumerics.h" ; then
       SiconosNumerics_inc="yes"
       if test "$ac_dir" != "/usr/include" ; then
	  SICONOSNUMERICS_INCLUDES="-I$ac_dir"
	  AC_MSG_RESULT([yes])
       fi       
       break	
    fi
done
if test "$SiconosNumerics_inc" = "" ; then
    AC_MSG_ERROR(Cannot find SICONOS Numerics header files (SiconosNumerics.h). Please install SICONOS Numerics.)
    exit 1;
fi
AC_SUBST(SICONOSNUMERICS_INCLUDES)

AC_MSG_CHECKING([SICONOS Numerics - libraries])
for ac_dir in		\
    /usr/lib		\
    /usr/local/lib	\
    $NUMERICSPATH/lib	\
    ; \
do  if test -r "$ac_dir/libNumericsSolvers_lcp.a" ; then
       SiconosNumerics_lib="yes"
       # SICONOSNUMERICS_LIBRARIES="-L$ac_dir -lNumericsSolvers_cfd -lNumericsSolvers_lcp -lNumericsSolvers_rp -lNumericsSolvers_cfp -lNumericsSolvers_rd -lodepack"
       SICONOSNUMERICS_LIBRARIES="-L$ac_dir -lNumericsSolvers_cfd -lNumericsSolvers_lcp -lNumericsSolvers_rp -lNumericsSolvers_cfp -lNumericsSolvers_rd"
       AC_MSG_RESULT([yes])
       break	
    fi
done
if test "$SiconosNumerics_lib" = "" ; then
    AC_MSG_ERROR(Cannot find SICONOS Numerics library (libNumericsSolvers_lcp.a). Please install SICONOS Numerics.)
    exit 1;
fi
AC_SUBST(SICONOSNUMERICS_LIBRARIES)




# Lapack: no include files (FORTRAN libraries)

AC_MSG_CHECKING([local lapack installation])
if test "$with_locallapack" = no -o "$with_locallapack" = yes -o "$with_locallapack" = ""; then
   AC_MSG_RESULT(option --with-locallapack not selected : installed lapack used)
   with_locallapack=no
else
   AC_MSG_RESULT(option  --with-locallapack selected :locally installed lapack used)
fi

AC_MSG_CHECKING([Lapack])

if test "$with_locallapack" = no; then
   list_dir="/usr/lib /usr/local/lib"
else
   list_dir="$with_locallapack $with_locallapack/lib"
fi

for ac_dir in $list_dir;

do  if test -r "$ac_dir/liblapack.so" ; then
       checkLibraryVersion $ac_dir/liblapack.so $LAPACK_VER
       if test $? -eq 0; then
	  AC_MSG_ERROR(Lapack library (liblapack.so) : the minimal required version is $LAPACK_VER)
	  exit 1;
       fi
       lapack_lib="yes"
       LAPACK_LIBRARIES="-L$ac_dir -llapack"
       AC_MSG_RESULT([yes])
       break	
    fi
done

# test static library
if test "$lapack_lib" = "" ; then
    for ac_dir in $list_dir;
    do  if test -r "$ac_dir/liblapack.a" ; then
	   lapack_lib="yes"
	   LAPACK_LIBRARIES="-L$ac_dir -llapack"
	   AC_MSG_RESULT([yes - found only a static version])
	   break
	else
	   AC_MSG_ERROR(Cannot find Lapack library. Please install Lapack.)
	   exit 1;
	fi
    done
fi
AC_SUBST(LAPACK_LIBRARIES)

# Lapack++

case "$target" in
    *-linux*)
      libname="liblapack++.so"
      ;;
    *-apple-darwin*)
      libname="liblapack++.dylib"  
      ;;
esac     

AC_MSG_CHECKING([local lapackpp installation])
if test "$with_locallapackpp" = no -o "$with_locallapackpp" = yes -o "$with_locallapackpp" = ""; then
   AC_MSG_RESULT(option --with-locallapackpp not selected : installed lapackpp used)
   with_locallapackpp=no
else
   AC_MSG_RESULT(option  --with-locallapackpp selected :locally installed lapackpp used)
fi

AC_MSG_CHECKING([Lapackpp - headers])

if test "$with_locallapackpp" = no; then
   list_dir="/usr/include /usr/local/include"
else
   list_dir="$with_locallapackpp/include"
fi

for ac_dir in $list_dir;

#for ac_dir in		\
#    /usr/include	\
#    /usr/local/include	\
#    ; \
do  if test -r "$ac_dir/lapackpp/lapackpp.h" ; then
       lapackpp_inc="yes"
       if test "$ac_dir" != "/usr/include" ; then
	  LAPACKPP_INCLUDES="-I$ac_dir/lapackpp"
	  AC_MSG_RESULT([yes])
       fi       
       break	
    fi
done
if test "$lapackpp_inc" = "" ; then
    AC_MSG_ERROR(Cannot find Lapack++ header files (lapackpp.h). Please install Lapack++.)
    exit 1;
fi
AC_SUBST(LAPACKPP_INCLUDES)

AC_MSG_CHECKING([Lapackpp - libraries])

if test "$with_locallapackpp" = no; then
   list_dir="/usr/lib /usr/local/lib"
else
   list_dir="$with_locallapackpp/lib"
fi

for ac_dir in $list_dir;
#for ac_dir in		\
#    /usr/lib		\
#    /usr/local/lib	\
#    ; \
do  if test -r "$ac_dir/$libname" ; then
       checkLibraryVersion $ac_dir/$libname $LAPACKPP_VER
       if test $? -eq 0; then
	  AC_MSG_ERROR(Lapack++ library ($libname) : the minimal required version is $LAPACKPP_VER)
	  exit 1;
       fi
       lapackpp_lib="yes"
       LAPACKPP_LIBRARIES="-L$ac_dir -llapack++ -llamatrix++ -lblas++"
       AC_MSG_RESULT([yes])
       break	
    fi
done
if test "$lapackpp_lib" = "" ; then
    AC_MSG_ERROR(Cannot find Lapack++ library ($libname). Please install Lapack++.)
    exit 1;
fi
AC_SUBST(LAPACKPP_LIBRARIES)

CXX_LAPACKPP_FLAG="-DLA_COMPLEX_SUPPORT" 
AC_SUBST(CXX_LAPACKPP_FLAG)

# Nana
AC_MSG_CHECKING([local nana installation])
if test "$with_localnana" = no -o "$with_localnana" = yes -o "$with_localnana" = ""; then
   AC_MSG_RESULT(option --with-localnana not selected : installed nana used)
   with_localnana=no
else
   AC_MSG_RESULT(option  --with-localnana selected :locally installed nana used)
fi

AC_MSG_CHECKING([nana - headers])

### TRY OF CONDITIONNAL CHECK
if test "$with_localnana" = no; then
   list_dir="/usr/include /usr/local/include"
else
   list_dir="$with_localnana/include"
fi
#echo list_dir is $list_dir
for ac_dir in $list_dir;
#for ac_dir in		\
#    /usr/include	\
#    /usr/local/include	\
#    ; \
do  if test -r "$ac_dir/nana.h" ; then
       nana_inc="yes"
       if test "$ac_dir" != "/usr/include" ; then
	  NANA_INCLUDES="-I$ac_dir"
	  AC_MSG_RESULT([yes])
       fi       
       break	
    fi
done
if test "$nana_inc" = "" ; then
    AC_MSG_ERROR(Cannot find nana header files (nana.h). Please install nana.)
    exit 1;
fi
AC_SUBST(NANA_INCLUDES)

AC_MSG_CHECKING([nana - libraries])
if test "$with_localnana" = no; then
   list_dir="/usr/lib /usr/local/lib"
else
   list_dir="$with_localnana/lib"
fi
for ac_dir in $list_dir;
#for ac_dir in		\
#    /usr/lib		\
#    /usr/local/lib	\
#    ; \
do  if test -r "$ac_dir/libnana.a" ; then
       nana_lib="yes"
       NANA_LIBRARIES="-L$ac_dir -lnana"
       AC_MSG_RESULT([yes])
       break	
    fi
done
if test "$nana_lib" = "" ; then
    AC_MSG_ERROR(Cannot find nana library (libnana.a). Please install nana.)
    exit 1;
fi
AC_SUBST(NANA_LIBRARIES)

CXX_WITHOUTNANA_FLAG="-DWITHOUT_NANA"

# Libxml2

AC_MSG_CHECKING([local libxml installation])
if test "$with_locallibxml" = no -o "$with_locallibxml" = yes -o "$with_locallibxml" = ""; then
   AC_MSG_RESULT(option --with-locallibxml not selected : installed libxml used)
   with_locallibxml=no
else
   AC_MSG_RESULT(option  --with-locallibxml selected :locally installed libxml used)
fi

AC_MSG_CHECKING([libxml2])

if test "$with_locallibxml" = no; then
   list_dir="/usr/bin /usr/local/bin"
else
   list_dir="$with_locallibxml/bin"
fi
echo list_dir is $list_dir
for ac_dir in $list_dir;

#for ac_dir in		\
#    /usr/bin		\
#    /usr/local/bin	\
#    ; \
do  if test -r "$ac_dir/xml2-config" ; then
       libxml2_config="yes"
# don't work on debussy, it uses the bad xml2-config which is in /usr/local !
#       if test "$ac_dir" != "/usr/bin" ; then
#	  LIBXML2_INCLUDES="`$ac_dir/xml2-config --cflags`"
#	  LIBXML2_LIBRARIES="`$ac_dir/xml2-config --libs`"
#       else 
#	  LIBXML2_INCLUDES="`xml2-config --cflags`" 
#	  LIBXML2_LIBRARIES="`xml2-config --libs`"
#       fi
       checkLibraryVersionNumber `$ac_dir/xml2-config --version` $LIBXML_VER
       if test $? -eq 0; then
	  AC_MSG_ERROR(Libxml library  : the minimal required version is $LIBXML_VER)
	  exit 1;
	  break;
       fi
       LIBXML2_INCLUDES="`$ac_dir/xml2-config --cflags`"
       LIBXML2_LIBRARIES="`$ac_dir/xml2-config --libs`"
       AC_MSG_RESULT([yes])       
       break	
    fi
done
#if test "$libxml2_config" = "" ; then
#    AC_MSG_ERROR(Cannot find libxml2 (xml2-config). Please install libxml2.)
#    exit 1;
#fi
# test static library
if test "$libxml2_config" = "" ; then
    if test "$with_locallibxml" = no; then
       list_dir="/usr/lib /usr/local/lib"
    else
       list_dir="$with_locallibxml/lib"
    fi
    echo list_dir is $list_dir
    for ac_dir in $list_dir;
    do  if test -r "$ac_dir/libxml2.a" ; then
	   libxml2_lib="yes"
	   LIBXML2_LIBRARIES="-L$ac_dir -llibxml2"
	   AC_MSG_RESULT([yes - found only a static version])
	   break
	else
	   AC_MSG_ERROR(Cannot find Libxml2 library. Please install Libxml2.)
	   exit 1;
	fi
    done
fi

AC_SUBST(LIBXML2_INCLUDES)
AC_SUBST(LIBXML2_LIBRARIES)

# cppunit

AC_MSG_CHECKING([if cppunit is enabled])
if test "$enable_cppunit" = no; then
   AC_MSG_RESULT(option --enable-cppunit not selected : Unable to run test suite)
else
   AC_MSG_RESULT(option --enable-cppunit selected : able to run test suite)
   
   AC_MSG_CHECKING([local cppunit installation])
   if test "$with_localcppunit" = no -o "$with_localcppunit" = yes -o "$with_localcppunit" = ""; then
      AC_MSG_RESULT(option --with-localcppunit not selected : installed cppunit used)
      with_localcppunit=no
   else
	AC_MSG_RESULT(option  --with-localcppunit selected :locally installed cppunit used)
   fi   

   AC_MSG_CHECKING([cppunit])
   
   if test "$with_localcppunit" = no; then
      list_dir="/usr/bin /usr/local/bin"
   else
	list_dir="$with_localcppunit/bin"
   fi
   #echo list_dir is $list_dir
   for ac_dir in $list_dir;

   #for ac_dir in		\
   #    /usr/bin		\
   #    /usr/local/bin	\
   #    ; \
   do  if test -r "$ac_dir/cppunit-config" ; then
	  checkLibraryVersionNumber `$ac_dir/cppunit-config --version` $CPPUNIT_VER
	  if test $? -eq 0; then
	     AC_MSG_ERROR(Cppunit library  : the minimal required version is $CPPUNIT_VER)
	     exit 1;
	  fi   
	  cppunit_config="yes"
	  CPPUNIT_INCLUDES="`$ac_dir/cppunit-config --cflags` -w"
	  CPPUNIT_LIBRARIES="`$ac_dir/cppunit-config --libs`"
	  AC_MSG_RESULT([yes])       
	  break	
       fi
   done
   if test "$cppunit_config" = "" ; then
      AC_MSG_ERROR(Cannot find cppunit (cppunit-config). Please install cppunit.)
      exit 1;
   fi
   AC_SUBST(CPPUNIT_INCLUDES)
   AC_SUBST(CPPUNIT_LIBRARIES)
fi



# COMPILE FLAGS

SICONOS_EXT_INCLUDES="$SICONOSNUMERICS_INCLUDES $NANA_INCLUDES $LAPACKPP_INCLUDES $BLAS_INCLUDES $LAPACK_INCLUDES $LIBXML2_INCLUDES"
AC_SUBST(SICONOS_EXT_INCLUDES)

case "$target" in
    *-linux*)
      SICONOS_EXT_LDFLAGS="$SICONOSNUMERICS_LIBRARIES $NANA_LIBRARIES $LAPACKPP_LIBRARIES $BLAS_LIBRARIES $LAPACK_LIBRARIES $LIBXML2_LIBRARIES"
      ;;
    *-apple-darwin*)
      SICONOS_EXT_LDFLAGS="$SICONOSNUMERICS_LIBRARIES $NANA_LIBRARIES $LAPACKPP_LIBRARIES $BLAS_LIBRARIES $LAPACK_LIBRARIES $LIBXML2_LIBRARIES -L/sw/lib -lg2c"
      ;;
esac
AC_SUBST(SICONOS_EXT_LDFLAGS)

CXXFLAGS_BASIC="-fbounds-check -w -Wno-deprecated -fPIC"
CXXFLAGS_DEBUG="-g"

case "$target" in
    *-linux*)
      CXXFLAGS="$CXXFLAGS_DEBUG $CXXFLAGS_BASIC $CXX_WITHOUTNANA_FLAG $CXX_LAPACKPP_FLAG"
      CXX_DYNLIB="-shared"
      ;;
    *-apple-darwin*)
      CXXFLAGS="$CXXFLAGS_DEBUG $CXXFLAGS_BASIC $CXX_WITHOUTNANA_FLAG $CXX_LAPACKPP_FLAG -flat_namespace"
      CXX_DYNLIB="-dynamiclib"
      ;;
esac
AC_SUBST(CXX_DYNLIB)

# Doxygen
DOXYGEN_CFG=$srcdir/config/DoxygenSiconos.cfg
AC_SUBST(DOXYGEN_CFG)

# output files
AC_CONFIG_FILES([Makefile
                 config/Makefile
		 python/Makefile
                 siconos/Makefile
                 src/Makefile
                 src/model/Makefile
                 src/model/test/Makefile
                 src/modelingTools/Makefile
                 src/simulationTools/Makefile
                 src/plugin/Makefile
                 src/utils/Makefile
		 src/utils/SiconosException/Makefile
                 src/utils/SiconosAlgebra/Makefile
                 src/utils/SiconosAlgebra/test/Makefile
                 src/utils/SiconosMemory/Makefile
                 src/utils/SiconosMemory/test/Makefile
                 src/utils/SiconosSharedLibrary/Makefile
                 src/xml/Makefile
                 src/xml/test/Makefile
		 sample/Makefile
		 sample/BouncingBall/Makefile
		 sample/RollingBalls/Makefile
		 sample/DoubleContact/Makefile
		 sample/UltraBalls/Makefile
		 sample/Ball2D/Makefile])
AC_OUTPUT
