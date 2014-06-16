#ifndef KERNELCONFIG_H
#define KERNELCONFIG_H
#define WITH_CMAKE

// Is cblas available? 
#cmakedefine HAS_CBLAS

// Where does it comes from? 
#cmakedefine HAS_MKL_CBLAS
#cmakedefine HAS_ACCELERATE // includes also lapack from Accelerate
#cmakedefine HAS_ATLAS_CBLAS 
#cmakedefine HAS_OpenBLAS // includes also lapacke from lapack/netlib

// Which Lapack? 
#cmakedefine HAS_MKL_LAPACKE
#cmakedefine HAS_ATLAS_LAPACK
#cmakedefine HAS_LAPACKE

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

// Say which version of C++ was used to compile the Kernel
#define KERNEL_CXXVERSION @CXX_VERSION@

#endif /*KERNELCONFIG_H*/

