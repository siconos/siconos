AC_DEFUN([ACX_CHECK_NUMERICS], [
AC_PREREQ(2.57)

# Numerics: no include files (FORTRAN libraries)
if test "$with_numerics" = no -o "$with_numerics" = yes -o "$with_numerics" = ""; then
    AC_MSG_RESULT(option --with-numerics not selected : default (/usr/local) numerics used)
    with_numerics=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local /sw /usr"  
      ;;
    *-linux*)
      list_dir="/usr /usr/local"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-numerics selected :locally installed numerics used)
   list_dir="$with_numerics"
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
	if test -r "$ac_dir/lib/libSiconosNumerics.$libsuffix" && test -r "$ac_dir/include/Siconos/SiconosNumerics.h" ; then
       		NUMERICS_INCLUDES="-I$ac_dir/include/Siconos"
       		NUMERICS_LIBRARIES="-L$ac_dir/lib -lSiconosNumerics"
       		NUMERICS_PATH="$ac_dir/"
       		numerics_lib="yes"
		dynlib="yes"
		echo "yes"
       		AC_MSG_RESULT([Library $ac_dir/lib/libSiconosNumerics.$libsuffix selected]) 
       		break;  
    	else
		echo "no"
	fi
done

# test static library
if test "$numerics_lib" = "no" ; then
    for ac_dir in $list_dir;
    do  AC_MSG_CHECKING([for libSiconosNumerics.a in $ac_dir])
        if test -r "$ac_dir/lib/libSiconosNumerics.a" ; then
	    	numerics_lib="yes"
		NUMERICS_INCLUDES="-I$ac_dir/include/Siconos"
       		NUMERICS_LIBRARIES="-L$ac_dir/lib -lSiconosNumerics"
       		NUMERICS_PATH="$ac_dir/"	    	   
		echo "yes"
	    	AC_MSG_RESULT([Library $ac_dir/lib/libNumerics.a selected])
	    	break;
    	else
		echo "no"
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


	     
])dnl ACX_CHECK_NUMERICS
