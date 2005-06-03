AC_DEFUN([ACX_CHECK_BLAS], [
AC_PREREQ(2.57)

# Blas: no include files (FORTRAN libraries)
if test "$with_localblas" = no -o "$with_localblas" = yes -o "$with_localblas" = ""; then
    AC_MSG_RESULT(option --with-localblas not selected : installed blas used)
    with_localblas=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local/lib /sw/lib /usr/lib "  
      ;;
    *-linux*)
      list_dir="/usr/lib /usr/local/lib"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-localblas selected :locally installed blas used)
   list_dir="$with_localblas $with_localblas/lib"
fi

case "$target" in
    *-apple-darwin*)
      libsuffix="dylib"
      ;;
    *-linux*)
      libsuffix="so"
      ;;
esac     


blas_lib="no"

for ac_dir in $list_dir;
do dynlib=no;    
   if test -r "$ac_dir/libblas.$libsuffix" ; then
       AC_MSG_CHECKING([for libblas.$libsuffix in $ac_dir])
       ACX_CHECK_DYNLIB_VER([$ac_dir/libblas.$libsuffix], [$BLAS_VER], [dynlib="yes"] ,[dynlib="higher version required"])	
       AC_MSG_RESULT($dynlib)
       if test "$dynlib" = "yes"; then
	   blas_lib="yes"
           BLAS_LIBRARIES="-L$ac_dir -lblas"
           break
       fi	
    fi
done

# test static library
if test "$blas_lib" = "no" ; then
    for ac_dir in $list_dir;
    do  AC_MSG_CHECKING([for libblas.a in $ac_dir])
        if test -r "$ac_dir/libblas.a" ; then
	    blas_lib="yes"
	    
	    BLAS_LIBRARIES="-L$ac_dir -lblas"
	    AC_MSG_RESULT([$blas_lib])
	    break
	fi
        AC_MSG_RESULT([$blas_lib])
    done
fi
# result of test
if test "$blas_lib" = "yes" ; then
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


	     
])dnl ACX_CHECK_BLAS
