AC_DEFUN([ACX_CHECK_CPPUNIT], [
AC_PREREQ(2.57)
# Libxml 2 
if test "$with_localcppunit" = no -o "$with_localcppunit" = yes -o "$with_localcppunit" = ""; then
    AC_MSG_RESULT(option --with-localcppunit not selected : installed cppunit used)
    with_locallibxml=no
   case "$target" in
    *-apple-darwin*)
      list_dir="/usr/local/bin /sw/bin /usr/bin "
      ;;
    *-linux*)
      list_dir="/usr/local/bin /usr/bin"  
      ;;
esac     
else
   AC_MSG_RESULT(option  --with-localcppunit selected :locally installed cppunit used)
   list_dir="$with_localcppunit/bin"
fi

cppunit_config="no"
for ac_dir in $list_dir;
do dynlib=no;  
   AC_MSG_CHECKING([for cppunit in $ac_dir])   
   if test -r "$ac_dir/cppunit-config" ; then
	echo Found cppunit $ac_dir/cppunit. version `$ac_dir/cppunit-config --version`
     	checkLibraryVersionNumber `$ac_dir/cppunit-config --version` $CPPUNIT_VER
       	if test $? -eq 0; then
	  AC_MSG_RESULT(but the minimal required version is $CPPUNIT_VER)
       	else         
          cppunit_config="yes"
	  CPPUNIT_INCLUDES="`$ac_dir/cppunit-config --cflags` -w"
	  CPPUNIT_LIBRARIES="`$ac_dir/cppunit-config --libs`"
          AC_MSG_RESULT([Good version - binary $ac_dir/cppunit selected]) 
	  break;  
       fi    	
    fi
done


if test "$cppunit_config" = "yes" ; then
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


	     
])dnl ACX_CHECK_CPPUNIT
