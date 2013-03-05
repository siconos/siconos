
// BLAS level 2
// symmetric matrices

//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <stddef.h>
#include <iostream>
#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;
namespace bindings = boost::numeric::bindings;
namespace tag = boost::numeric::bindings::tag;

using std::cout;
using std::cin;
using std::endl; 

typedef double real_t; 
typedef ublas::vector<real_t> vct_t;

typedef ublas::symmetric_matrix<
  real_t, ublas::upper, ublas::column_major
> ucsymm_t; 
typedef ublas::symmetric_matrix<
  real_t, ublas::lower, ublas::column_major
> lcsymm_t; 
typedef ublas::symmetric_matrix<
  real_t, ublas::upper, ublas::row_major
> ursymm_t; 
typedef ublas::symmetric_matrix<
  real_t, ublas::lower, ublas::row_major
> lrsymm_t; 

typedef ublas::matrix<real_t, ublas::column_major> cm_t;
typedef ublas::matrix<real_t, ublas::row_major> rm_t;
typedef ublas::symmetric_adaptor<cm_t, ublas::upper> ucsa_t; 
typedef ublas::symmetric_adaptor<cm_t, ublas::lower> lcsa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::upper> ursa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::lower> lrsa_t; 

int main (int argc, char **argv) {
  size_t n = 0;
  if (argc > 1) {
    n = atoi(argv [1]);
  }
  
  cout << endl; 
  if (n <= 0) {
    cout << "n -> ";
    cin >> n;
  }
  cout << endl; 

  vct_t vx (n), vy (n); 
  blas::set (1., vx); 
  print_v (vx, "vx"); 
  vy (1) = 1.; 
  print_v (vy, "vy"); 

  // symmetric matrix

  ucsymm_t ucs (n, n); 
  lcsymm_t lcs (n, n); 
  ursymm_t urs (n, n); 
  lrsymm_t lrs (n, n); 

  init_symm (ucs, 'u'); 
  init_symm (lcs, 'l'); 
  init_symm (urs, 'u'); 
  init_symm (lrs, 'l'); 

  print_m (ucs, "ucs");
  cout << endl; 
  print_m (lcs, "lcs");
  cout << endl; 
  print_m (urs, "urs");
  cout << endl; 
  print_m (lrs, "lrs");
  cout << endl; 

  // m += x x^T 
  blas::spr (1.0, vx, ucs); 
  print_m (ucs, "ucs += x x^T"); 
  cout << endl; 
  blas::spr (1.0, vx, lcs); 
  print_m (lcs, "lcs += x x^T"); 
  cout << endl; 
  blas::spr (1.0, vx, urs); 
  print_m (urs, "urs += x x^T"); 
  cout << endl; 
  blas::spr (1.0, vx, lrs); 
  print_m (lrs, "lrs += x x^T"); 
  cout << endl; 

  // m += x y^T + y x^T
  blas::spr2 (1.0, vx, vy, ucs); 
  print_m (ucs, "ucs += x y^T + y x^T"); 
  cout << endl; 
  blas::spr2 (1.0, vx, vy, lcs); 
  print_m (lcs, "lcs += x y^T + y x^T"); 
  cout << endl; 
  blas::spr2 (1.0, vx, vy, urs); 
  print_m (urs, "urs += x y^T + y x^T"); 
  cout << endl; 
  blas::spr2 (1.0, vx, vy, lrs); 
  print_m (lrs, "lrs += x y^T + y x^T"); 
  cout << endl; 


  // symmetric adaptor

  cm_t cmu (n, n); 
  cm_t cml (n, n); 
  rm_t rmu (n, n); 
  rm_t rml (n, n); 

  ucsa_t ucsa (cmu); 
  lcsa_t lcsa (cml); 
  ursa_t ursa (rmu); 
  lrsa_t lrsa (rml); 

  init_symm (ucsa, 'u'); 
  init_symm (lcsa, 'l'); 
  init_symm (ursa, 'u'); 
  init_symm (lrsa, 'l'); 

  print_m (ucsa, "ucsa");
  cout << endl; 
  print_m (lcsa, "lcsa");
  cout << endl; 
  print_m (ursa, "ursa");
  cout << endl; 
  print_m (lrsa, "lrsa");
  cout << endl; 

  // m += x x^T 
  blas::syr (1.0, vx, ucsa); 
  print_m (ucsa, "ucsa += x x^T"); 
  cout << endl; 
  blas::syr (1.0, vx, lcsa); 
  print_m (lcsa, "lcsa += x x^T"); 
  cout << endl; 
  blas::syr (1.0, vx, ursa); 
  print_m (ursa, "ursa += x x^T"); 
  cout << endl; 
  blas::syr (1.0, vx, lrsa); 
  print_m (lrsa, "lrsa += x x^T"); 
  cout << endl; 

  // m += x y^T + y x^T
  blas::syr2 (1.0, vx, vy, ucsa); 
  print_m (ucsa, "ucsa += x y^T + y x^T"); 
  cout << endl; 
  blas::syr2 (1.0, vx, vy, lcsa); 
  print_m (lcsa, "lcsa += x y^T + y x^T"); 
  cout << endl; 
  blas::syr2 (1.0, vx, vy, ursa); 
  print_m (ursa, "ursa += x y^T + y x^T"); 
  cout << endl; 
  blas::syr2 (1.0, vx, vy, lrsa); 
  print_m (lrsa, "lrsa += x y^T + y x^T"); 
  cout << endl; 

}
