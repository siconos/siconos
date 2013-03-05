
// solving A * X = B
// using driver function gesv()

#include <cstddef>
#include <iostream>
#include <vector>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#if !defined(TEST_MATLIB_UBLAS) && !defined(TEST_MATLIB_GLAS) && !defined(TEST_MATLIB_MTL) && !defined(TEST_MATLIB_EIGEN)
#define TEST_MATLIB_UBLAS
#endif

#if defined(TEST_MATLIB_UBLAS)

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
namespace ublas = boost::numeric::ublas;
typedef ublas::matrix<double, ublas::column_major> m_t;
typedef std::size_t size_type;

#elif defined(TEST_MATLIB_GLAS)

#include <boost/numeric/bindings/glas/dense_matrix.hpp>
#include <glas/toolbox/la/algorithm/operators.hpp>
using namespace glas::la;
typedef glas::dense_matrix<double, glas::column_orientation> m_t;
typedef std::ptrdiff_t size_type;

#elif defined(TEST_MATLIB_MTL)

#include <boost/numeric/bindings/mtl/dense2D.hpp>
#include <boost/numeric/mtl/operation/operators.hpp>
typedef mtl::dense2D<double, mtl::matrix::parameters<mtl::tag::col_major> > m_t;
typedef std::ptrdiff_t size_type;

#elif defined(TEST_MATLIB_EIGEN)

#include <boost/numeric/bindings/eigen/matrix.hpp>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_t;
typedef int size_type;

#endif

#include "utils.h"
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;
using std::cout;
using std::endl;

int main() {

  cout << endl; 

  size_type n = 5;   
  m_t a (n, n);   // system matrix 

  size_type nrhs = 2; 
  m_t x (n, nrhs), b (n, nrhs);  // b -- right-hand side matrix

  init_symm (a); 
  //     [n   n-1 n-2  ... 1]
  //     [n-1 n   n-1  ... 2]
  // a = [n-2 n-1 n    ... 3]
  //     [        ...       ]
  //     [1   2   ...  n-1 n]

  m_t aa (a); // copy of a, because a is `lost' after gesv()

#if defined(TEST_MATLIB_UBLAS)
  ublas::matrix_column<m_t> xc0 (x, 0), xc1 (x, 1); 
  for (int i = 0; i < xc0.size(); ++i) {
    xc0 (i) = 1.;
    xc1 (i) = 2.; 
  }
  b = prod (a, x); 
#elif defined(TEST_MATLIB_GLAS) || defined(TEST_MATLIB_MTL) || defined(TEST_MATLIB_EIGEN)
  for (int i = 0; i < bindings::size_row (x); ++i) {
    x (i,0) = 1.;
    x (i,1) = 2.; 
  }
  b = a * x;
#endif

  print_m (a, "A"); 
  cout << endl; 
  print_m (b, "B"); 
  cout << endl; 

//  lapack::gesv (a, b);  // solving the system, b contains x 
//  no ipiv less version is currently provided, so fall back to using ipiv
  std::vector<fortran_int_t> ipiv(n);
  lapack::gesv (a, ipiv, b);  // solving the system, b contains x 

  print_m (b, "X");
  cout << endl; 

#if defined(TEST_MATLIB_UBLAS)
  x = prod (aa, b); 
#elif defined(TEST_MATLIB_GLAS) || defined(TEST_MATLIB_MTL) || defined(TEST_MATLIB_EIGEN)
  x = aa * b;
#endif
  print_m (x, "B = A X"); 

  cout << endl; 

}

