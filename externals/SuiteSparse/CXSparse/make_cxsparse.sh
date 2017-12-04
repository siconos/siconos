#!/bin/sh
set -e

if ! [ -e SuiteSparse-5.0.0.tar.gz ]; then
    curl -O -L 'http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.0.0.tar.gz'
fi
tar -xzf SuiteSparse-5.0.0.tar.gz
make -C SuiteSparse/SuiteSparse_config

cat preamble.txt >cxsparse_di.c
cat $(ls SuiteSparse/CXSparse/Source/*.c | grep -v cs_convert.c) >>cxsparse_di.c
echo cxsparse_di.c created.
cat preamble.txt >cxsparse_dl.c
echo '#define CS_LONG' >>cxsparse_dl.c
echo '#include "cxsparse_di.c"' >>cxsparse_dl.c
echo cxsparse_dl.c created.
cat preamble.txt >cxsparse_ci.c
#cat preamble_convert.txt >>cxsparse_ci.c
echo '#include "SiconosConfig.h"' >>cxsparse_ci.c
echo '#ifndef BUILD_AS_CPP' >>cxsparse_ci.c
echo '#define CS_COMPLEX' >>cxsparse_ci.c
echo '#include "cxsparse_di.c"' >>cxsparse_ci.c
echo '#endif' >>cxsparse_ci.c
echo cxsparse_ci.c created.
cat preamble.txt >cxsparse_cl.c
#cat preamble_convert.txt >>cxsparse_cl.c
echo '#include "SiconosConfig.h"' >>cxsparse_cl.c
echo '#ifndef BUILD_AS_CPP' >>cxsparse_cl.c
echo '#define CS_COMPLEX' >>cxsparse_cl.c
echo '#define CS_LONG' >>cxsparse_cl.c
echo '#include "cxsparse_di.c"' >>cxsparse_cl.c
echo '#endif' >>cxsparse_cl.c
echo cxsparse_cl.c created.
#cat preamble.txt preamble_convert.txt >cxsparse_conv.c
cat preamble.txt >cxsparse_conv.c
echo '#include "SiconosConfig.h"' >>cxsparse_conv.c
echo '#ifndef BUILD_AS_CPP' >>cxsparse_conv.c
cat SuiteSparse/CXSparse/Source/cs_convert.c >>cxsparse_conv.c
echo '#endif' >>cxsparse_conv.c
echo cxsparse_convert.c created.

cp -v SuiteSparse/CXSparse/Include/cs.h .
cp -v SuiteSparse/include/SuiteSparse_config.h .
echo cs.h and SuiteSparse_config.h copied.

cp -v SuiteSparse/LICENSE.txt .
