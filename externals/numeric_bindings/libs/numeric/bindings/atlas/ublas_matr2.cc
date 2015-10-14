
// BLAS level 2

//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <iostream>
#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include "utils.h" 

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;

using std::cout;
using std::endl; 

typedef ublas::vector<double> vct_t;
typedef ublas::matrix<double, ublas::row_major> rm_t;
typedef ublas::matrix<double, ublas::column_major> cm_t;

int main() {

  cout << endl; 

  vct_t vx (2);
  vct_t vy (4); 

  // row major matrix
  rm_t rm (4, 2);
  init_m (rm, kpp (1)); 
  print_m (rm, "row major matrix m"); 
  cout << endl; 

  blas::set (1., vx);
  print_v (vx, "vx"); 
  cout << endl; 

  // vy = m vx
  blas::gemv ( 1.0, rm, vx, 0.0, vy);
  print_v (vy, "vy = m vx"); 
  cout << endl; 

  blas::set (1., vy); 
  print_v (vy, "vy"); 
  cout << endl; 

  // vx = m^T vy
  blas::gemv ( 1.0, bindings::trans(rm), vy, 0.0, vx);
  print_v (vx, "vx = m^T vy"); 
  cout << endl; 

  cout << endl; 

  // column major matrix
  cm_t cm (4, 2);
  init_m (cm, kpp (1)); 
  print_m (cm, "column major matrix m"); 
  cout << endl; 

  blas::set (1., vx);
  print_v (vx, "vx"); 
  cout << endl; 

  // vy = m vx
  blas::gemv (1.0, cm, vx, 0.0, vy);
  print_v (vy, "vy = m vx"); 
  cout << endl; 

  blas::set (1., vy); 
  print_v (vy, "vy"); 
  cout << endl; 

  // vx = m^T vy
  blas::gemv ( 1.0, bindings::trans(cm), vy, 0.0, vx);
  print_v (vx, "vx = m^T vy"); 
  cout << endl; 

}
