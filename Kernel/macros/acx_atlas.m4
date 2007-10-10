# CHECK_ATLAS

AC_DEFUN([ACX_ATLAS],[

AC_ARG_ENABLE([atlas],
	AS_HELP_STRING([--enable-atlas],[enable ATLAS (default=yes)]),
  	enable_atlas="$enableval",
	enable_atlas="yes")

if test "$enable_atlas" = "no"; then
    return
fi


AC_ARG_WITH([atlas-include-dir],
	AS_HELP_STRING([atlas-include-dir=DIR],
		[specify the directory for atlas headers (optional)]),
		[if test x"$withval" != x""; then
	        	ac_atlas_include_dir="$withval"
		 fi])

AC_ARG_WITH([atlas-lib-dir],
	AS_HELP_STRING([atlas-lib-dir=DIR],
		[specify the directory for atlas libraries (optional)]),
		[if test x"$withval" != x""; then
	        	ac_atlas_lib_dir="$withval"
		 fi])


dnl first we check the system location for atlas libraries
if test -n "$ac_atlas_include_dir"; then
	ATLAS_CPPFLAGS="-I$ac_atlas_include_dir"
fi

if test -n "$ac_atlas_lib_dir"; then
	ATLAS_LDFLAGS="-L$ac_atlas_lib_dir"
fi

if test -z $ATLAS_CPPFLAGS; then
	for ac_atlas_path_tmp in /usr /usr/local /opt ; do
		if test -d "$ac_atlas_path_tmp/include/atlas" && test -r "$ac_atlas_path_tmp/include/atlas"; then			
		    ATLAS_CPPFLAGS="-I$ac_atlas_path_tmp/include/atlas"
		    break;	
		fi
	done
fi

if test -z $ATLAS_LDFLAGS; then
	for ac_atlas_path_tmp in /usr /usr/local /opt ; do
		if test -d "$ac_atlas_path_tmp/lib/atlas" && test -r "$ac_atlas_path_tmp/lib/atlas"; then			
		    ATLAS_LDFLAGS="-L$ac_atlas_path_tmp/lib/atlas"
		    break;	
		fi
	done
fi

CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$ATLAS_CPPFLAGS $CPPFLAGS"
export CPPFLAGS

LDFLAGS_SAVED="$LDFLAGS"
LDFLAGS="$ATLAS_LDFLAGS $LDFLAGS"
export LDFLAGS

acx_atlas_ok=yes
      
AC_CHECK_HEADERS([atlas_enum.h],,acx_atlas_ok=no )
AC_CHECK_HEADERS([clapack.h],,acx_atlas_ok=no)
AC_CHECK_HEADERS([cblas.h],,acx_atlas_ok=no)

AC_CHECK_LIB([atlas], [ATL_buildinfo], , [acx_atlas_ok=no])

AC_SEARCH_LIBS([clapack_dgetrf],[lapack_atlas lapack],,[acx_atlas_ok=no])

AC_SUBST(ATLAS_CPPFLAGS)
AC_SUBST(ATLAS_LDFLAGS)

#LDFLAGS is needed now for other AC_* as -llapack_atlas or -latlas
#may be in some non-std dir (i.e /usr/lib/atlas)
#LDFLAGS=$LDFLAGS_SAVED CPPFLAGS=$CPPFLAGS_SAVED

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_atlas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_ATLAS,1,[Define if you have a ATLAS library.]),[$1])
        :
else
        acx_atlas_ok=no
        $2
fi

])dnl ACX_CHECK_ATLAS

