#!/bin/sh
set -e

# if ! [ -e SuiteSparse-5.0.0.tar.gz ]; then
#     curl -O -L 'http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-5.0.0.tar.gz'
# fi
# tar -xzf SuiteSparse-5.0.0.tar.gz
# make -C SuiteSparse/SuiteSparse_config


#----------- CXSparse
mkdir -p CXSparse
cat preamble.txt > CXSparse/cxsparse_di.c
cat $(ls SuiteSparse/CXSparse/Source/*.c | grep -v cs_convert.c) >> CXSparse/cxsparse_di.c
echo CXSparse/cxsparse_di.c created.
cat preamble.txt >CXSparse/xsparse_dl.c
echo '#define CS_LONG' >>CXSparse/cxsparse_dl.c
echo '#include "cxsparse_di.c"' >>CXSparse/cxsparse_dl.c
echo CXSparse/cxsparse_dl.c created.
cat preamble.txt >CXSparse/cxsparse_ci.c
#cat preamble_convert.txt >>CXSparse/cxsparse_ci.c
echo '#include "SiconosConfig.h"' >>CXSparse/cxsparse_ci.c
echo '#ifndef BUILD_AS_CPP' >>CXSparse/cxsparse_ci.c
echo '#define CS_COMPLEX' >>CXSparse/cxsparse_ci.c
echo '#include "cxsparse_di.c"' >>CXSparse/cxsparse_ci.c
echo '#endif' >>CXSparse/cxsparse_ci.c
echo CXSparse/cxsparse_ci.c created.
cat preamble.txt >CXSparse/cxsparse_cl.c
#cat preamble_convert.txt >>CXSparse/cxsparse_cl.c
echo '#include "SiconosConfig.h"' >>CXSparse/cxsparse_cl.c
echo '#ifndef BUILD_AS_CPP' >>CXSparse/cxsparse_cl.c
echo '#define CS_COMPLEX' >>CXSparse/cxsparse_cl.c
echo '#define CS_LONG' >>CXSparse/cxsparse_cl.c
echo '#include "cxsparse_di.c"' >>CXSparse/cxsparse_cl.c
echo '#endif' >>CXSparse/cxsparse_cl.c
echo CXSparse/cxsparse_cl.c created.
#cat preamble.txt preamble_convert.txt >CXSparse/cxsparse_conv.c
cat preamble.txt >CXSparse/cxsparse_conv.c
echo '#include "SiconosConfig.h"' >>CXSparse/cxsparse_conv.c
echo '#ifndef BUILD_AS_CPP' >>CXSparse/cxsparse_conv.c
cat SuiteSparse/CXSparse/Source/cs_convert.c >>CXSparse/cxsparse_conv.c
echo '#endif' >>CXSparse/cxsparse_conv.c
echo CXSparse/cxsparse_convert.c created.

cp -v SuiteSparse/CXSparse/Include/cs.h CXSparse/
echo cs.h copied.

#----------- LDL
mkdir -p LDL

cat preamble.txt >LDL/ldl_di.c
echo '#define LDL_LONG' >> LDL/ldl_di.c
cat $(ls SuiteSparse/LDL/Source/*.c | grep -v cs_convert.c) >>LDL/ldl_di.c
echo LDL/ldl_di.c created.


cp -v SuiteSparse/LDL/Include/ldl.h LDL
echo ldl.h copied.

cp -v SuiteSparse/include/SuiteSparse_config.h .


cd CXSparse
ln -sf ../SuiteSparse_config.h SuiteSparse_config.h
cd ../LDL
ln -sf ../SuiteSparse_config.h SuiteSparse_config.h
cd ..

echo SuiteSparse_config.h copied.

cp -v SuiteSparse/LICENSE.txt .
