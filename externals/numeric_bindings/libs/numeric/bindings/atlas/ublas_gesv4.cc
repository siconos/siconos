
// solving A * X = B
// using driver function gesv()
// with ublas::vector<> as RHS

#include <cstddef>
#include <iostream>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include <boost/numeric/bindings/std/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp> 

namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;

using std::size_t; 
using std::cout;
using std::endl; 

typedef ublas::matrix<double, ublas::column_major> m_t;
typedef ublas::vector<double> v_t; 

int main() {

  cout << endl; 
  size_t n = 3;

  m_t a (n, n);   // system matrix 
  a(0,0) = 1.; a(0,1) = 1.; a(0,2) = 1.;
  a(1,0) = 2.; a(1,1) = 3.; a(1,2) = 1.;
  a(2,0) = 1.; a(2,1) = -1.; a(2,2) = -1.;

  v_t b (n);  // right-hand side vector
  b(0) = 4.; b(1) = 9.; b(2) = -2.; 

  cout << "A: " << a << endl; 
  cout << "B: " << b << endl; 

  std::vector< int > pivot( bindings::size1( a ) );
  lapack::gesv (a, pivot, b);  
  cout << "X: " << b << endl; 
  cout << endl; 

}

