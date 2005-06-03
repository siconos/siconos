AC_DEFUN([ACX_CHECK_LAPACK], [
AC_PREREQ(2.57)

# Lapack: no include files (FORTRAN libraries)
if test "$with_locallapack" = no -o "$with_locallapack" = yes -o "$with_locallapack" = ""; then
    AC_MSG_RESULT(option --with-locallapack not selected : installed lapack used)
    with_locallapack=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local/lib /sw/lib /usr/lib "  
      ;;
    *-linux*)
      list_dir="/usr/lib /usr/local/lib"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-locallapack selected :locally installed lapack used)
   list_dir="$with_locallapack $with_locallapack/lib"
fi

case "$target" in
    *-apple-darwin*)
      libsuffix="dylib"
      ;;
    *-linux*)
      libsuffix="so"
      ;;
esac  


lapack_lib="no"

for ac_dir in $list_dir;
do dynlib=no;    
   if test -r "$ac_dir/liblapack.$libsuffix" ; then
       AC_MSG_CHECKING([for liblapack.$libsuffix in $ac_dir])
       ACX_CHECK_DYNLIB_VER([$ac_dir/liblapack.$libsuffix], [$LAPACK_VER], [dynlib="yes"] ,[dynlib="higher version required"])	
       AC_MSG_RESULT($dynlib)
       if test "$dynlib" = "yes"; then
	   lapack_lib="yes"
           LAPACK_LIBRARIES="-L$ac_dir -llapack"
           break
       fi	
    fi
done

# test static library
if test "$lapack_lib" = "no" ; then
    for ac_dir in $list_dir;
    do  AC_MSG_CHECKING([for liblapack.a in $ac_dir])
        if test -r "$ac_dir/liblapack.a" ; then
	    lapack_lib="yes"
	    
	    LAPACK_LIBRARIES="-L$ac_dir -llapack"
	    AC_MSG_RESULT([$lapack_lib])
	    break
	fi
        AC_MSG_RESULT([$lapack_lib])
    done
fi
# result of test
if test "$lapack_lib" = "yes" ; then
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
