/*
 * 
 * Copyright (c) Kresimir Fresl 2002 
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Author acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_CBLAS_H
#define BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_CBLAS_H



#ifdef HAS_OpenBLAS
// OpenBlas use a different signature for some functions.
// See for example https://bitbucket.org/blaze-lib/blaze/issues/50/cannot-convert-const-std-complex-to-const.
#define OPENBLAS_CONST_FLOAT_CAST(x) reinterpret_cast<const float*>(x)
#define OPENBLAS_CONST_DOUBLE_CAST(x) reinterpret_cast<const double*>(x)
#define OPENBLAS_FLOAT_CAST(x) reinterpret_cast<float*>(x)
#define OPENBLAS_DOUBLE_CAST(x) reinterpret_cast<double*>(x)
#define OPENBLAS_OPENBLAS_COMPLEX_FLOAT_CAST(x) reinterpret_cast<openblas_complex_float*>(x)
#define OPENBLAS_OPENBLAS_COMPLEX_DOUBLE_CAST(x) reinterpret_cast<openblas_complex_double*>(x)

#else
#define OPENBLAS_CONST_FLOAT_CAST(x) x
#define OPENBLAS_CONST_DOUBLE_CAST(x) x
#define OPENBLAS_FLOAT_CAST(x) x
#define OPENBLAS_DOUBLE_CAST(x) x
#define OPENBLAS_OPENBLAS_COMPLEX_FLOAT_CAST(x) x
#define OPENBLAS_OPENBLAS_COMPLEX_DOUBLE_CAST(x) x
#endif

//
// MKL-specific CBLAS include
//
#if defined BOOST_NUMERIC_BINDINGS_BLAS_MKL

extern "C" {
#include <mkl_cblas.h>
//#include <mkl_service.h>
//
// mkl_types.h defines P4 macro which breaks MPL, undefine it here.
//
#undef P4
}

//
// Default CBLAS include
//
#else

extern "C" {
#include <cblas.h>
}

#endif
#endif 
