AC_DEFUN([ACX_CHECK_LAPACKPP], [
AC_PREREQ(2.57)

# Lapackpp: no include files (FORTRAN libraries)
if test "$with_locallapackpp" = no -o "$with_locallapackpp" = yes -o "$with_locallapackpp" = ""; then
    AC_MSG_RESULT(option --with-locallapackpp not selected : installed lapackpp used)
    with_locallapackpp=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local /sw /usr"  
      ;;
    *-linux*)
      list_dir="/usr /usr/local"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-locallapackpp selected :locally installed lapackpp used)
   list_dir="$with_locallapack $with_locallapackpp"
fi

case "$target" in
    *-apple-darwin*)
      libsuffix="dylib"
      ;;
    *-linux*)
      libsuffix="so"
      ;;
esac  


lapackpp_lib="no"

for ac_dir in $list_dir;
do dynlib=no;    
   if test -r "$ac_dir/lib/liblapackpp.$libsuffix" && test -r "$ac_dir/include/lapackpp/lapackpp.h" ; then
       AC_MSG_CHECKING([for liblapackpp.$libsuffix in $ac_dir])
       checkLibraryVersion $ac_dir/lib/liblapackpp.$libsuffix $LAPACKPP_VER
	if test $? -eq 0; then
	  AC_MSG_RESULT(but the minimal required version is $LAPACKPP_VER)
       else         
          LAPACKPP_INCLUDES="-I$ac_dir/include/lapackpp"
          LAPACKPP_LIBRARIES="-L$ac_dir -llapack++ -llamatrix++ -lblas++"
          lapackpp_lib="yes"
          AC_MSG_RESULT([Good version - library $ac_dir/lib/lapackpp.$libsuffix selected]) 
	  break;  
       fi 
    fi
done

# test static library
#if test "$lapack_lib" = "no" ; then
#    for ac_dir in $list_dir;
#    do  AC_MSG_CHECKING([for liblapack.a in $ac_dir])
#        if test -r "$ac_dir/liblapack.a" ; then
#	    lapack_lib="yes"
	    
#	    LAPACK_LIBRARIES="-L$ac_dir -llapack"
#	    AC_MSG_RESULT([$lapack_lib])
#	    break
#	fi
#        AC_MSG_RESULT([$lapack_lib])
#    done
#fi
# result of test
if test "$lapackpp_lib" = "yes" ; then
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
