#! /bin/sh

sed "s/SWIG_append_msg/SWIG_append_errmsg/g" pySiconos_wrap.cxx > pySiconos_wrap.cxx.res ;
mv pySiconos_wrap.cxx.res pySiconos_wrap.cxx; 
exit 0;
