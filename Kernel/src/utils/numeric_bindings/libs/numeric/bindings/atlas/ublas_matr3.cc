
// BLAS level 3

//#define F_USE_STD_VECTOR

//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <iostream>
#include <boost/numeric/bindings/blas/level3.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/trans.hpp>
#ifdef F_USE_STD_VECTOR
#include <vector>
#include <boost/numeric/bindings/std/vector.hpp> 
#endif 
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;

using std::cout;
using std::endl; 

#ifndef F_USE_STD_VECTOR
typedef ublas::matrix<double, ublas::row_major> m_t;
#else
typedef ublas::matrix<double, ublas::column_major, std::vector<double> > m_t;
#endif 

int main() {

  cout << endl; 

  m_t a (4, 4);
  init_m (a, kpp (1)); 
  print_m (a, "a"); 
  cout << endl; 

  m_t b (4, 6);
  init_m (b, cls1()); 
  print_m (b, "b"); 
  cout << endl; 
  
  m_t c (4, 6);

  // c = a b
  blas::gemm ( 1.0, a, b, 0.0, c); 
  print_m (c, "c = a b"); 
  cout << endl; 
  blas::gemm ( 1.0, a, b, 0.0, c); 
  print_m (c, "c = a b"); 
  cout << endl; 

  init_m (c, const_val<double> (1)); 
  print_m (c, "c"); 
  cout << endl; 
  // c = 2 a b + 0.5 c
  blas::gemm (2.0, a, b, 0.05, c);
  print_m (c, "c = 2 a b + 0.05 c"); 
  cout << endl; 

  m_t d (6, 4);

  // d = b^T a^T
  blas::gemm ( 1.0, bindings::trans(b), bindings::trans(a), 0.0, d);
  print_m (d, "d = b^T a^T"); 
  cout << endl; 

  // c = a^T b 
  blas::gemm (1.0, bindings::trans(a), b, 0.0, c); 
  print_m (c, "c = a^T b"); 
  cout << endl; 

  // d = b^T a
  blas::gemm ( 1.0, bindings::trans(b), a, 0.0, d);
  print_m (d, "d = b^T a"); 
  cout << endl; 

  init_m (d, const_val<double> (0)); 
  ublas::matrix_range<m_t> br (b, ublas::range (0, 4), ublas::range (0, 4)); 
  ublas::matrix_range<m_t> dr (d, ublas::range (1, 5), ublas::range (0, 4)); 

  // d[1..5][0..4] = a b[0..4][0..4]  
  blas::gemm ( 1.0, a, br, 0.0, dr); 
  print_m (d, "d[1..5][0..4] = a b[0..4][0..4]"); 
  cout << endl; 
  
  // d[1..5][0..4] = b[0..4][0..4] a
  blas::gemm ( 1.0, br, a, 0.0, dr); 
  print_m (d, "d[1..5][0..4] = b[0..4][0..4] a"); 
  cout << endl; 
  
  // d[1..5][0..4] = b[0..4][0..4] a^T
  blas::gemm ( 1.0, br, bindings::trans(a), 0.0, dr);
  print_m (d, "d[1..5][0..4] = b[0..4][0..4] a^T"); 
  cout << endl; 

  // d[1..5][0..4] = a b[0..4][0..4]^T
  blas::gemm ( 1.0, a, bindings::trans(br), 0.0, dr);
  print_m (d, "d[1..5][0..4] = a b[0..4][0..4]^T"); 
  cout << endl; 

}
