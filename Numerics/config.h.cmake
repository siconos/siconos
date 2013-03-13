#ifndef NUMERICSCONFIG_H
#define NUMERICSCONFIG_H
#define WITH_CMAKE

#cmakedefine HAVE_PATHFERRIS
#cmakedefine HAVE_MLCPSIMPLEX
#cmakedefine HAVE_TIME_H
#cmakedefine HAVE_SYSTIMES_H
#cmakedefine WITH_MUMPS
#cmakedefine WITH_TIMERS
#cmakedefine DUMP_PROBLEM
#cmakedefine WITH_FCLIB
#cmakedefine BUILD_AS_CPP


// Is cblas available? 
#cmakedefine HAS_CBLAS

// Where does it comes from? 
#cmakedefine HAS_MKL_CBLAS
#cmakedefine HAS_ACCELERATE // includes also lapack from Accelerate
#cmakedefine HAS_ATLAS_CBLAS 
#cmakedefine HAS_OpenBLAS // it *MAY* also includes lapacke or lapack from netlib
#cmakedefine HAS_GenericCBLAS

// Which Lapack? 
#cmakedefine HAS_MKL_LAPACKE
#cmakedefine HAS_ATLAS_LAPACK
#cmakedefine HAS_LAPACKE

// Which functions are defined in lapack? 
#cmakedefine HAS_LAPACK_DGESVD
#cmakedefine HAS_LAPACK_DTRTRS
#cmakedefine HAS_LAPACK_DGELS
#endif
