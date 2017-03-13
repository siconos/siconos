#ifndef SICONOSCONFIG_H
#define SICONOSCONFIG_H

#define WITH_CMAKE
#cmakedefine HAVE_SICONOS_KERNEL
#cmakedefine HAVE_SICONOS_IO
#cmakedefine HAVE_SICONOS_MECHANICS
#cmakedefine HAVE_SICONOS_CONTROL
#cmakedefine HAVE_PATHFERRIS
#cmakedefine HAVE_PATHVI
#cmakedefine HAVE_MLCPSIMPLEX
#cmakedefine HAVE_TIME_H
#cmakedefine HAVE_SYSTIMES_H
#cmakedefine HAVE_MPI
#cmakedefine WITH_MUMPS
#cmakedefine WITH_UMFPACK
#cmakedefine WITH_SUPERLU
#cmakedefine WITH_SUPERLU_MT
#cmakedefine WITH_SUPERLU_dist
#cmakedefine WITH_MKL_PARDISO
#cmakedefine WITH_MKL_SPBLAS
#cmakedefine SICONOS_MKL_32
#cmakedefine SICONOS_MKL_64
#cmakedefine WITH_TIMERS
#cmakedefine DUMP_PROBLEM
#cmakedefine WITH_FCLIB
#cmakedefine BUILD_AS_CPP
#cmakedefine WITH_LPSOLVE
#cmakedefine HAS_EXTREME_POINT_ALGO
#cmakedefine HAVE_VTK
#cmakedefine HAVE_BULLET
#cmakedefine HAVE_OCC
#cmakedefine HAVE_SERIALIZATION
#cmakedefine WITH_SERIALIZATION
#cmakedefine HAVE_GENERATION
#cmakedefine WITH_HDF5
#cmakedefine WITH_OPENMP

// Is cblas available? 
#cmakedefine HAS_CBLAS

// Where does cblas comes from? 
#cmakedefine HAS_MKL_CBLAS
#cmakedefine HAS_ACCELERATE // includes also lapack from Accelerate
#cmakedefine HAS_ATLAS_CBLAS 
#cmakedefine HAS_OpenBLAS // it *MAY* also includes lapacke or lapack from netlib
#cmakedefine HAS_GenericCBLAS

// Which Lapack? 
#cmakedefine HAS_MKL_LAPACKE
#cmakedefine HAS_ATLAS_LAPACK
#cmakedefine HAS_MATLAB_LAPACK
#cmakedefine HAS_LAPACKE // lapacke.h has been found
#cmakedefine HAS_CLAPACK  // clapack.h has been found
#cmakedefine HAS_OpenBLAS_LAPACK

// Which functions are defined in lapack? 
#cmakedefine HAS_LAPACK_DGESVD
#cmakedefine HAS_LAPACK_DTRTRS
#cmakedefine HAS_LAPACK_DGELS

// Some definitions required for boost numeric_bindings
#if defined(HAS_CBLAS)
#define BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
#endif

#if defined(HAS_MKL_CBLAS)
#define BOOST_NUMERIC_BINDINGS_BLAS_MKL
#endif

// Gams stuff
#cmakedefine GAMS_MODELS_SOURCE_DIR "@GAMS_MODELS_SOURCE_DIR@"
#cmakedefine GAMS_MODELS_SHARE_DIR "@GAMS_MODELS_SHARE_DIR@"
#cmakedefine GAMS_DIR "@GAMS_DIR@"
#cmakedefine HAVE_GAMS_C_API

// Which version of C++ was used to compile siconos, needed for swig
#define SICONOS_CXXVERSION @CXXVERSION@
#cmakedefine SICONOS_USE_BOOST_FOR_CXX11
#cmakedefine SICONOS_USE_MAP_FOR_HASH
#cmakedefine SICONOS_STD_SHARED_PTR
#cmakedefine SICONOS_STD_ARRAY
#cmakedefine SICONOS_STD_UNORDERED_MAP
#cmakedefine SICONOS_STD_TUPLE
#cmakedefine SICONOS_STD_TO_STRING
#cmakedefine SICONOS_STD_FUNCTIONAL

// are int 64 bits longs
#cmakedefine SICONOS_INT64

// use to force 32 bits int when creating numpy array
// Useful to support old scipy version (< 0.14.0)
#cmakedefine SICONOS_FORCE_NPY_INT32
#endif
