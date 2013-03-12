* Which blas implementations are Siconos-Numerics compliant?
  * Netlib 
  * Apple framework accelerate
  * OpenBlas
  * MKL intel
  * Atlas 
  
* Which blas implementations are Siconos-Numerics compliant?
  * Netlib  (lapacke)
  * Apple framework accelerate
  * OpenBlas (with lapacke)
  * MKL intel (lapacke again)
  * Atlas but since it's not complete, it's probably better to use lapacke from netlib, with atlas-cblas.

* About Cblas :
  
  Blas implementation is based on cblas. All interfaces are the same. Routines are cblas_XXX.
  For a proper documentation with examples, see for instance : 
  
  - http://software.intel.com/sites/products/documentation/hpc/mkl/mklman/index.htm#GUID-A491FDF5-ACA0-4EDC-8F74-00306D0EBE05.htm
  - https://developer.apple.com/library/mac/#documentation/Accelerate/Reference/BLAS_Ref/Reference/reference.html

  For Siconos Numerics everything is defined in SiconosBlas.h
  
  ** Libraries and headers:
  The searched libraries depends on the blas implementation. We used the FindBlas.cmake from cmake with a few changes (for la)

* About Lapack
  It seems that lapacke interface is the common "point" to several implementations. 
  See :http://www.netlib.org/lapack/lapacke.html
  And http://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/

  Apple and atlas are a bit different and based on (old?) clapack. 
  
  Everything is defined in SiconosLapack.h with implementation-specific devel in SiconosAtlas.h, SiconosApple.h and SiconosLapacke.h

* Which libraries and headers are searched by cmake? 

  Openblas : libopenblas.so .dylib ..., lapacke.h
  mkl : lib depends on intel version, mkl_cblas.h, mkl_lapacke.h
  netlib : libcblas, cblas.h, lapacke.h
  atlas : lib? , cblas.h, clapack.h
  apple : Framework accelerate, Accelerate.h, clapack.h, cblas.h
  
  

* Difference between lapack implem
  * char */int*  in Apple, char in lapacke, atlas
  * some routines are not implemented in the basic atlas (dgesvd, dtrtrs, dgels)
  * char are (sometimes) replaced by enum in atlas
  * no more work vectors in high-level lapacke routines (dgels ...) while they are still present in apple or clapack.h 

These difference are handled with WRAP_XXX macros in SiconosLapack.h.


    
