.. _siconos_dependencies:

Siconos required and optional dependencies
==========================================

* a compiler suite, with c++, c and gfortran compilers.

* cmake (version > 2.8.7, 3.x will be better) (https://cmake.org)

* boost (http://www.boost.org)

* blas and lapack (see :ref:`about_blas_lapack`)

  
To generate the documentation, you will need :

* doxygen
* sphinx

For the python bindings:

* python (>= 2.7)
* swig (>= 2.0)

To run tests:

* cppunit


.. _about_blas_lapack:

About blas and Lapack
=====================

The BLAS (Basic Linear Algebra Subprograms, http://www.netlib.org/blas/) are routines that provide standard building blocks for performing basic vector and matrix operations, while LAPACK (http://www.netlib.org/lapack/#_presentation) provides routines for solving systems of simultaneous linear equations, least-squares solutions of linear systems of equations, eigenvalue problems, and singular value problems.
Different implementations are available, such as:

* openblas (http://www.openblas.net),
* the one from MKL (https://software.intel.com/en-us/intel-mkl),
* Accelerate framework on Macosx (https://developer.apple.com/library/prerelease/mac/documentation/Accelerate/Reference/BLAS_Ref/index.html) ...
  
For siconos we recommand:

* accelerate on Macosx
* OpenBLAS + lapacke on linux systems

Warning : we do not provide support for atlas.

Anyway, power users can still choose their favorite BLAS and LAPACK(E) vendor.

Blas, lapack setup of your system will be checked during cmake call.

If the process failed or if you need a specific implementation, the following variables may be provided to cmake to help the searching process (see :ref:`siconos_detailed_install`)

* BLA_VENDOR : Blas implementation type (used also as hint for lapack).
  One of :
  * OpenBLAS
  * Matlab
  * Intel10_32 (intel mkl v10 32 bit)
  * Intel10_64lp (intel mkl v10+ 64 bit, threaded code, lp64 model)
  * Intel10_64lp_seq (intel mkl v10+ 64 bit, sequential code, lp64 model)
  * Intel10_64ilp (intel mkl v10+ 64 bit, threaded code, ilp64 model)
  * Intel10_64ilp_seq (intel mkl v10+ 64 bit, sequential code, ilp64 model)
  * Apple
  * Generic

* BLAS_ROOT : Blas implementation location.
* LAPACK_ROOT : Lapack implementation location.


About Boost
===========

Boost provides a lot of useful C++ binaries, especially Ublas, a C++ template class library that provides BLAS level 1, 2, 3 functionalities 
for dense, packed and sparse matrices.

Ublas is used in Siconos for matrices and vectors definition and implementation.

About Boost: http://www.boost.org/

About Ublas: http://www.boost.org/libs/numeric/ublas/doc/index.htm

Install (note that an adequate Boost version comes with most linux distributions and thus no more install is required.)

To know how to get and install Boost, see 
Boost Getting Started.

Note that we also use boost-bindings:
"Boost Bindings is a bindings library (not just) for Boost.Ublas. It offers an easy way of calling BLAS, LAPACK, UMFPACK, MUMPS and many other mature legacy numerical codes from within C++."

They are distributed and installed with the Siconos but you can also get the last version here: 
http://mathema.tician.de/software/boost-bindings

GMP
===

"GMP is a free library for arbitrary precision arithmetic, operating on signed integers, rational numbers, and floating point numbers ... "

This library usually comes with gcc. If not see http://gmplib.org/ for download and installation instructions.
