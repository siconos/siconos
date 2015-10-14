/* ------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.     */
/* Davis.  All Rights Reserved.  See ../README for License.                  */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.           */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                      */
/* ------------------------------------------------------------------------- */


/***********************************************************************/
/*         UMFPACK Copyright, License and Availability                 */
/***********************************************************************/
/*
 *
 * UMFPACK Version 4.1 (Apr. 30, 2003),  Copyright (c) 2003 by Timothy A.
 * Davis.  All Rights Reserved.
 *
 * UMFPACK License:
 *
 *   Your use or distribution of UMFPACK or any modified version of
 *   UMFPACK implies that you agree to this License.
 *
 *   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 *   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 *
 *   Permission is hereby granted to use or copy this program, provided
 *   that the Copyright, this License, and the Availability of the original
 *   version is retained on all copies.  User documentation of any code that
 *   uses UMFPACK or any modified version of UMFPACK code must cite the
 *   Copyright, this License, the Availability note, and "Used by permission."
 *   Permission to modify the code and to distribute modified code is granted,
 *   provided the Copyright, this License, and the Availability note are
 *   retained, and a notice that the code was modified is included.  This
 *   software was developed with support from the National Science Foundation,
 *   and is provided to you free of charge.
 *
 * Availability:
 *
 *   http://www.cise.ufl.edu/research/sparse/umfpack
 *
 */

/* Used by permission. */


/* Simple demo program for UMFPACK               */
/* from UMFPACK Version 4.1 Quick Start Guide    */

/* modified by Kresimir Fresl, 2003              */
/* UMFPACK bindings & ublas::compressed_matrix<> */


#if !defined(TEST_MATLIB_UBLAS) && !defined(TEST_MATLIB_GLAS) && !defined(TEST_MATLIB_MTL) && !defined(TEST_MATLIB_EIGEN)
#define TEST_MATLIB_UBLAS
#endif

#include <iostream>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#if defined(TEST_MATLIB_UBLAS)

#include <boost/numeric/bindings/ublas/matrix_sparse.hpp>
namespace ublas = boost::numeric::ublas;
typedef ublas::compressed_matrix<double, ublas::column_major, 0, 
    ublas::unbounded_array<int>, ublas::unbounded_array<double> > m_t;
typedef m_t& mi_t;
#define AA(i,j, val) a(i,j) = val

#elif defined(TEST_MATLIB_GLAS)

#include <boost/numeric/bindings/glas/compressed.hpp>
typedef glas::sparse_matrix< double, glas::compressed_sparse_structure<
    glas::column_orientation, std::ptrdiff_t, std::ptrdiff_t, 0> > m_t;
typedef m_t& mi_t;
#define AA(i,j, val) insert(a,i,j, val)

#elif defined(TEST_MATLIB_MTL)

#include <boost/numeric/bindings/mtl/compressed2D.hpp>
typedef mtl::compressed2D<double, mtl::matrix::parameters<mtl::tag::col_major> > m_t;
typedef mtl::matrix::inserter<m_t> mi_t;
#define AA(i,j, val) a(i,j) = val

#elif defined(TEST_MATLIB_EIGEN)

#include <boost/numeric/bindings/eigen/sparsematrix.hpp>
typedef Eigen::SparseMatrix<double> m_t;
typedef Eigen::RandomSetter<m_t> mi_t;
#define AA(i,j, val) a(i,j) = val

#endif

typedef boost::numeric::ublas::vector<double> v_t;
namespace umf = boost::numeric::bindings::umfpack;

int main() {

  m_t A (5,5
#if !defined(TEST_MATLIB_EIGEN)
        ,12
#endif
        );
  v_t B (5), X (5);

  {
    mi_t a(A);
    AA(0,0, 2.); AA(0,1, 3.);
    AA(1,0, 3.); AA(1,2, 4.); AA(1,4, 6.);
    AA(2,1,-1.); AA(2,2,-3.); AA(2,3, 2.);
    AA(3,2, 1.);
    AA(4,1, 4.); AA(4,2, 2.); AA(4,4, 1.);
  }
  B(0) = 8.; B(1) = 45.; B(2) = -3.; B(3) = 3.; B(4) = 19.;

  umf::symbolic_type<double> Symbolic;
  umf::numeric_type<double> Numeric;

  umf::symbolic (A, Symbolic);
  umf::numeric (A, Symbolic, Numeric);
  umf::solve (A, X, B, Numeric);

  std::cout << X << std::endl;  // output: [5](1,2,3,4,5)
}
