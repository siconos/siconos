#! /bin/sh

LIBSICONOS="$SICONOSPATH/src/utils/SiconosException/libSiconosException.a $SICONOSPATH/src/utils/SiconosSharedLibrary/libSiconosSharedLibrary.a $SICONOSPATH/src/utils/NewSiconosVector/libSiconosMatrixVector.a $SICONOSPATH/src/utils/SiconosMemory/libSiconosMemory.a $SICONOSPATH/src/xml/libSiconosXml.a $SICONOSPATH/src/modelformalisation/libModelFormalisation.a $SICONOSPATH/src/modelstrategy/libModelStrategy.a $SICONOSPATH/src/model/libModel.a"

#LIBNUMERICS="$NUMERICSPATH/lib/libNumericsSolvers_cfd.a $NUMERICSPATH/lib/libNumericsSolvers_cfp.a $NUMERICSPATH/lib/libNumericsSolvers_lcp.a $NUMERICSPATH/lib/libNumericsSolvers_rd.a $NUMERICSPATH/lib/libNumericsSolvers_rp.a $NUMERICSPATH/lib/libodepack.a"

LIBNUMERICS="$NUMERICSPATH/lib/libNumericsSolvers_cfd.a $NUMERICSPATH/lib/libNumericsSolvers_cfp.a $NUMERICSPATH/lib/libNumericsSolvers_lcp.a $NUMERICSPATH/lib/libNumericsSolvers_rd.a $NUMERICSPATH/lib/libNumericsSolvers_rp.a"

TARGET=libSiconos

mkdir foo
cd foo
for lib in ${LIBSICONOS} ${LIBNUMERICS}; do
    ar x $lib
done
cd ..
