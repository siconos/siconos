AC_DEFUN([ACX_CHECK_NUMERICS], [
AC_PREREQ(2.57)

# Numerics: no include files (FORTRAN libraries)
if test "$with_localnumerics" = no -o "$with_localnumerics" = yes -o "$with_localnumerics" = ""; then
    AC_MSG_RESULT(option --with-localnumerics not selected : installed numerics used)
    with_localnumerics=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local /sw /usr"  
      ;;
    *-linux*)
      list_dir="/usr /usr/local"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-localnumerics selected :locally installed numerics used)
   list_dir="$with_locallapack $with_localnumerics"
fi

case "$target" in
    *-apple-darwin*)
      libsuffix="dylib"
      ;;
    *-linux*)
      libsuffix="so"
      ;;
esac  


numerics_lib="no"
dynlib=no;    
for ac_dir in $list_dir;
do  AC_MSG_CHECKING([for libSiconosNumerics.$libsuffix in $ac_dir])
	   if test -r "$ac_dir/lib/libSiconosNumerics.$libsuffix" && test -r "$ac_dir/include/SiconosNumerics.h" ; then
       		NUMERICS_INCLUDES="-I$ac_dir/include/"
       		NUMERICS_LIBRARIES="-L$ac_dir/lib -lSiconosNumerics"
       		numerics_lib="yes"
		dynlib="yes"
       		AC_MSG_RESULT([Library $ac_dir/lib/libSiconosNumerics.$libsuffix selected]) 
       	break;  
    fi
done

# test static library
if test "$numerics_lib" = "no" ; then
    for ac_dir in $list_dir;
    do  AC_MSG_CHECKING([for libSiconosNumerics.a in $ac_dir])
        if test -r "$ac_dir/lib/libSiconosNumerics.a" ; then
	    	numerics_lib="yes"
		NUMERICS_INCLUDES="-I$ac_dir/include/"
       		NUMERICS_LIBRARIES="-L$ac_dir/lib -lSiconosNumerics"	    	   
	    	AC_MSG_RESULT([Library $ac_dir/lib/libNumerics.a selected])
	    break
	fi
    done
fi

# result of test
if test "$numerics_lib" = "yes" ; then
    if test "$dynlib" = "no"; then
       result="static version found"
       $1
    else
       result="dynamic version found" 
       $1
    fi                
else
    result="no"
    $2
fi


	     
])dnl ACX_CHECK_LAPACK
