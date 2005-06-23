AC_DEFUN([ACX_CHECK_NANA], [
AC_PREREQ(2.57)
echo
# Nana 2 
if test "$with_localnana" = no -o "$with_localnana" = yes -o "$with_localnana" = ""; then
    AC_MSG_RESULT(option --with-localnana not selected : installed nana used)
    with_localnana=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local /sw /usr "
      ;;
    *-linux*)
      list_dir="/usr/local /usr"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-localnana selected :locally installed nana used)
   list_dir="$with_localnana/"
fi

case "$target" in
    *-apple-darwin*)
      libsuffix="dylib"
      ;;
    *-linux*)
      libsuffix="so"
      ;;
esac  

nana_lib="no"
for ac_dir in $list_dir;
do dynlib=no;  
   AC_MSG_CHECKING([for libnana.$libsuffix in $ac_dir/lib])   
   if test -r "$ac_dir/lib/libnana.$libsuffix"  && test -r "$ac_dir/include/nana.h" ;then
          NANA_INCLUDES="$ac_dir/lib/nana.h"
          NANA_LIBRARIES="$ac_dir/lib/libnana.$libsuffix"
          nana_lib="yes"
          dynlib = "yes" 
          AC_MSG_RESULT([Library $ac_dir/lib/nana.$libsuffix selected]) 
	  break;  
   fi    	
done

# test static library
if test "$nana_lib" = "no" ; then
    for ac_dir in $list_dir;
    do  AC_MSG_CHECKING([for libnana.a in $ac_dir/lib])
        if test -r "$ac_dir/lib/libnana.a"  && test -r "$ac_dir/include/nana.h" ;  then
	    nana_lib="yes"
            NANA_INCLUDES="-I$ac_dir/include"       
	    NANA_LIBRARIES="-L$ac_dir/lib -lnana"
	    AC_MSG_RESULT([Found  a static library $ac_dir/lib/libnana.a])
            
        break
	fi
    done
fi



# result of test
if test "$nana_lib" = "yes" ; then
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


	     
])dnl ACX_CHECK_NANA
