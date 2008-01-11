# TRY TO FIND ATLAS INSTALLATION LAYOUT

AC_DEFUN([ACX_ATLAS],[

AC_ARG_ENABLE([atlas],
	AS_HELP_STRING([--enable-atlas],[enable ATLAS (default=yes)]),
  	enable_atlas="$enableval",
	enable_atlas="yes")

if test "$enable_atlas" = "no"; then
    return
fi

AC_ARG_WITH(atlas-include,,
	[  --with-atlas-include=<dir>   Use this ATLAS include directory],
        [case $withval in 
           yes | "") ;;
           *) ATLAS_CPPFLAGS="$ATLAS_CPPFLAGS -I$withval/atlas -I$withval";; 
         esac])

AC_ARG_WITH(atlas-lib,,
	[  --with-atlas-lib=<dir>   Use this ATLAS library directory],
        [case $withval in 
           yes | "") ;;
           *) ATLAS_LDFLAGS="$ATLAS_LDFLAGS -L$withval";;
         esac])  


# first we check the system location for atlas libraries
if test -n "$ac_atlas_include_dir"; then
	ATLAS_CPPFLAGS="$ATLAS_CPPFLAGS -I$ac_atlas_include_dir/atlas -I$ac_atlas_include_dir"
fi

if test -n "$ac_atlas_lib_dir"; then
	ATLAS_LDFLAGS="$ATLAS_LDFLAGS -L$ac_atlas_lib_dir"
fi

if test -z $ATLAS_CPPFLAGS; then
	for ac_atlas_path_tmp in /usr /usr/local /opt ; do
# as found with some pkg (deb, rpm)
		if test -d "$ac_atlas_path_tmp/include/atlas" && test -r "$ac_atlas_path_tmp/include/atlas"; then
		    if test -z $ATLAS_LDFLAGS; then 
			ATLAS_LDFLAGS="-L$ac_atlas_path_tmp/lib/atlas"
		    fi
		    ATLAS_CPPFLAGS="-I$ac_atlas_path_tmp/include/atlas"
		    break;
# as found with a default installation layout (atlas 3.8)
                else if test -d "$ac_atlas_path_tmp/atlas/include/atlas" && test -r "$ac_atlas_path_tmp/atlas/include/atlas"; then
                    if test -z $ATLAS_LDFLAGS; then 
			ATLAS_LDFLAGS="-L$ac_atlas_path_tmp/atlas/lib/atlas -L$ac_atlas_path_tmp/atlas/lib"
		    fi
		    ATLAS_CPPFLAGS="-I$ac_atlas_path_tmp/atlas/include/atlas -I$ac_atlas_path_tmp/atlas/include"
		    break; 	
		  fi
                fi
	done
fi

CPPFLAGS_SAVED="$CPPFLAGS"
LDFLAGS_SAVED="$LDFLAGS"
LIBS_SAVED="$LIBS"
FLIBS_SAVED="$FLIBS"

CPPFLAGS="$CPPFLAGS $ATLAS_CPPFLAGS"
export CPPFLAGS


LDFLAGS="$LDFLAGS $ATLAS_LDFLAGS"
export LDFLAGS

LIBS=
FLIBS=

acx_atlas_ok=yes

# header check
# debian with deb atlas 3.6 / install from atlas 3.8 => none of *.h have the same names...
# so we cannot run AC_CHECK_HEADERS on them

# necessary but not sufficient
AC_CHECK_HEADERS([clapack.h],,acx_atlas_ok=no)

AC_CHECK_HEADERS([cblas.h],,acx_atlas_ok=no)

# lib check
AC_CHECK_LIB([atlas], [ATL_buildinfo], , [acx_atlas_ok=no])

AC_F77_FUNC(fatlas_caxpby)
AC_CHECK_LIB([f77blas], [$fatlas_caxpby])

AC_CHECK_LIB([cblas], [catlas_daxpby], , [AC_CHECK_LIB([blas], [catlas_daxpy])])

AC_CHECK_LIB([lapack], [ATL_dgetrf])

AC_CHECK_LIB([lapack_atlas], [ATL_dgetri])

AC_SUBST(ATLAS_CPPFLAGS)
AC_SUBST(ATLAS_LDFLAGS)

ATLAS_INCLUDES="$ATLAS_CPPFLAGS"
ATLAS_LIBS="$LIBS"

# come back to defaults
LIBS="$LIBS_SAVED"
FLIBS="$FLIBS_SAVED"
LDFLAGS="$LDFLAGS_SAVED"
CPPFLAGS="$CPPFLAGS_SAVED"

# finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_atlas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_ATLAS,1,[Define if you have a ATLAS library.]),[$1])
        :
else
        acx_atlas_ok=no
        $2
fi
])dnl ACX_ATLAS

