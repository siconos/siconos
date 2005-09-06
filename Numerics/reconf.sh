#! /bin/sh
rm -f config.cache
libtoolize
aclocal -I ./m4 > /dev/null 
autoconf
autoheader
automake --add-missing
exit 0
