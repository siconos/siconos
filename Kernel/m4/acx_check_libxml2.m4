AC_DEFUN([ACX_CHECK_LIBXML2], [
AC_PREREQ(2.57)
echo
# Libxml 2 
if test "$with_locallibxml" = no -o "$with_locallibxml" = yes -o "$with_locallibxml" = ""; then
    AC_MSG_RESULT(option --with-locallibxml not selected : installed libxml used)
    with_locallibxml=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local /sw /usr "
      ;;
    *-linux*)
      list_dir="/usr/local /usr"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-locallibxml selected :locally installed libxml used)
   list_dir="$with_locallibxml/"
fi

case "$target" in
    *-apple-darwin*)
      libsuffix="dylib"
      ;;
    *-linux*)
      libsuffix="so"
      ;;
esac  

libxml_lib="no"


for ac_dir in $list_dir;
do dynlib=no;  
   AC_MSG_CHECKING([for libxml2.$libsuffix in $ac_dir/lib])   
   if test -r "$ac_dir/bin/xml2-config" ; then
	echo Found a dynamic library $ac_dir/lib/libxml2.$libsuffix. version `$ac_dir/bin/xml2-config --version`
        checkLibraryVersionNumber `$ac_dir/bin/xml2-config --version` $LIBXML_VER
       if test $? -eq 0; then
	  AC_MSG_RESULT(but the minimal required version is $LIBXML_VER)
       else         
          LIBXML2_INCLUDES="`$ac_dir/bin/xml2-config --cflags`"
          LIBXML2_LIBRARIES="`$ac_dir/bin/xml2-config --libs`"
	  if test "${LOCAL_LIB_PATH}test" = "test"; then
	    export LOCAL_LIB_PATH="$ac_dir/lib";
          else
	    export LOCAL_LIB_PATH=$LOCAL_LIB_PATH:"$ac_dir/lib";
	  fi;
          libxml_lib="yes"
          AC_MSG_RESULT([Good version - library $ac_dir/lib/libxml2.$libsuffix selected]) 
	  break;  
       fi    	
    fi
done



# test static library
#if test "$libxml_lib" = "no" ; then
#    for ac_dir in $list_dir;
#    do  AC_MSG_CHECKING([for libxml2.a in $ac_dir/lib])
#        if test -r "$ac_dir/lib/libxml2.a" ; then
#	    libxml_lib="yes"          
#	    LIBXML2_LIBRARIES="-L$ac_dir/lib -llibxml2"
#        break
#	fi
#	AC_MSG_RESULT([yes - Found  a static library $ac_dir/lib/$libxml_lib])
#    done
#fi
# result of test
if test "$libxml_lib" = "yes" ; then
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


	     
])dnl ACX_CHECK_LIBXML2
