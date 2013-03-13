#ifndef KERNELCONFIG_H
#define KERNELCONFIG_H
#define WITH_CMAKE
#define XML_SCHEMA "@CMAKE_INSTALL_PREFIX@/share/@PROJECT_PACKAGE_NAME@/SiconosModelSchema-V1.2.xsd"

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

#endif /*KERNELCONFIG_H*/

