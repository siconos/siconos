AC_DEFUN([ACX_CHECK_KERNEL], [
AC_PREREQ(2.57)

# Kernel
if test "$with_localkernel" = no -o "$with_localkernel" = yes -o "$with_localkernel" = ""; then
    AC_MSG_RESULT(option --with-localkernel not selected : installed kernel used)
    with_localkernel=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local /sw /usr"  
      ;;
    *-linux*)
      list_dir="/usr /usr/local"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-localkernel selected :locally installed kernel used)
   list_dir="$with_localkernel"
fi

case "$target" in
    *-apple-darwin*)
      libsuffix="dylib"
      ;;
    *-linux*)
      libsuffix="so"
      ;;
esac  


kernel_lib="no"
dynlib=no;    
for ac_dir in $list_dir;
do  AC_MSG_CHECKING([for libSiconosKernel.$libsuffix in $ac_dir]) 
	   if test -r "$ac_dir/lib/libSiconosKernel.$libsuffix" && test -r "$ac_dir/include/KernelDefaultConfig.h" ; then
       		KERNEL_INCLUDES="-I$ac_dir/include/"
       		KERNEL_LIBRARIES="-L$ac_dir/lib -lSiconosKernel"
       		KERNEL_PATH="$ac_dir/"
       		kernel_lib="yes"
		dynlib="yes"
		echo "yes"
       		AC_MSG_RESULT([Library $ac_dir/lib/libSiconosKernel.$libsuffix selected]) 
       	break;  
	else
	echo "no"
    fi
done

# test static library
if test "$kernel_lib" = "no" ; then
    for ac_dir in $list_dir;
    do  AC_MSG_CHECKING([for libSiconosKernel.a in $ac_dir])
        if test -r "$ac_dir/lib/libSiconosKernel.a" ; then
	    	kernel_lib="yes"
		KERNEL_INCLUDES="-I$ac_dir/include/"
       		KERNEL_LIBRARIES="-L$ac_dir/lib -lSiconosKernel"
       		KERNEL_PATH="$ac_dir/"	  
		echo "yes"  	   
	    	AC_MSG_RESULT([Library $ac_dir/lib/libSiconosKernel.a selected])
	    break
	else
	echo "no"
	fi
    done
fi

# result of test
if test "$kernel_lib" = "yes" ; then
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
