# CHECK_ATLAS

AC_DEFUN([ACX_CHECK_ATLAS],[

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
		    if test -z $ATLAS_LDFLAGS; then 
			ATLAS_LDFLAGS="-L$ac_atlas_path_tmp/lib/atlas"
		    fi
		    ATLAS_CPPFLAGS="-I$ac_atlas_path_tmp/include/atlas"
		    break;	
		fi
	done
fi

CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $ATLAS_CPPFLAGS"
export CPPFLAGS

LDFLAGS_SAVED="$LDFLAGS"
LDFLAGS="$LDFLAGS $ATLAS_LDFLAGS"
export LDFLAGS

acx_atlas_ok=yes
      
AC_CHECK_HEADERS([atlas_enum.h],,acx_atlas_ok=no )
AC_CHECK_HEADERS([clapack.h],,acx_atlas_ok=no)
AC_CHECK_HEADERS([cblas.h],,acx_atlas_ok=no)

AC_CHECK_LIB([atlas], [ATL_buildinfo], , [acx_atlas_ok=no])

AC_SUBST(ATLAS_CPPFLAGS)
AC_SUBST(ATLAS_LDFLAGS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_atlas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_ATLAS,1,[Define if you have a ATLAS library.]),[$1])
        :
else
        acx_atlas_ok=no
        $2
fi

])dnl ACX_CHECK_ATLAS

