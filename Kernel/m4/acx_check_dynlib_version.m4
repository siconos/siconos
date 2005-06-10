dnl @synopsis ACX_CHECK_DYNLIB_VER([LIB-PATH, VERSION-NUMBER, ACTION-IF-FOUND, ACTION-IF-NOT-FOUND])

AC_DEFUN([ACX_CHECK_DYNLIB_VER], [
AC_PREREQ(2.57)

# function checkLibraryVersion ()
# $1 : absolute path of the library (ex. : /usr/lib/liblapack.so)
# $2 : minimum version (ex. : 3.0.3)

lib=$1; 
ver=$2; 
lslib=`ls -l $lib`
vers=`echo $lslib | awk '{ print $NF}' |  awk 'BEGIN {FS = "."; } { printf "%d", ($[3] * 1000 + $[4]) * 1000 + $[5];}'`
verReq=`expr $ver  | awk 'BEGIN {FS = "."; } { printf "%d", ($[1] * 1000 + $[2]) * 1000 + $[3];}'`
if test -n "$vers" && test "$vers" -ge $verReq ; then
     $3;
else
     $4;
fi	        
])dnl ACX_CHECK_DYNLIB_VER

dnl @synopsis ACX_CHECK_VER([VERSION-NUMBER, VERSION-NUMBER-REQUIRED, ACTION-IF-FOUND, ACTION-IF-NOT-FOUND])

AC_DEFUN([ACX_CHECK_VER], [
AC_PREREQ(2.57)
# function checkLibraryVersionNumber ()
# $1 : version number ex. : 3.2.1)
# $2 : minimum version (ex. : 3.0.3)

vers=`expr $1 |  awk 'BEGIN {FS = "."; } { printf "%d", ($[1] * 1000 + $[2]) * 1000 + $[3];}'`
verReq=`expr $2  | awk 'BEGIN {FS = "."; } { printf "%d", ($[1] * 1000 + $[2]) * 1000 + $[3];}'`
if test -n "$vers" && test "$vers" -ge $verReq ; then
    $3;
else
    $4;
fi	     
])dnl ACX_CHECK_DYNLIB_VER


