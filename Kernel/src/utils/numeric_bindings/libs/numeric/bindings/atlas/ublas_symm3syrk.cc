
// BLAS level 3
// symmetric matrices, syrk 

#include <stddef.h>
#include <iostream>
#include <boost/numeric/bindings/blas/level3.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/upper.hpp>
#include <boost/numeric/bindings/lower.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;
namespace bindings = boost::numeric::bindings;

using std::cout;
using std::cin;
using std::endl; 

typedef double real_t; 

typedef ublas::matrix<real_t, ublas::column_major> cm_t;
typedef ublas::matrix<real_t, ublas::row_major> rm_t;
typedef ublas::symmetric_adaptor<cm_t, ublas::upper> ucsa_t; 
typedef ublas::symmetric_adaptor<cm_t, ublas::lower> lcsa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::upper> ursa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::lower> lrsa_t; 

int main (int argc, char **argv) {
  int n = 0, k = 0;
  if (argc > 1) {
    n = atoi(argv [1]);
  }
  if (argc > 2) {
    k = atoi(argv [2]);
  }

  if (n <= 0) {
    cout << "n -> ";
    cin >> n;
  }
  if (k <= 0) {
    cout << "k -> ";
    cin >> k; 
  }

  cm_t ac (n, k); 
  rm_t ar (n, k); 
  init_m (ac, rws1()); 
  init_m (ar, rws1());
  print_m (ac, "ac"); 
  cout << endl; 
  print_m (ar, "ar"); 
  cout << endl << endl;

  cm_t cmu (n, n); 
  cm_t cml (n, n); 
  rm_t rmu (n, n); 
  rm_t rml (n, n); 
  ucsa_t ucsa (cmu); 
  lcsa_t lcsa (cml); 
  ursa_t ursa (rmu); 
  lrsa_t lrsa (rml); 

  blas::syrk (1.0, ac, 0.0, ucsa); 
  blas::syrk (1.0, ac, 0.0, lcsa); 
  blas::syrk (1.0, ar, 0.0, ursa); 
  blas::syrk (1.0, ar, 0.0, lrsa); 

  print_m (ucsa, "ucsa");
  cout << endl; 
  print_m (lcsa, "lcsa");
  cout << endl; 
  print_m (ursa, "ursa");
  cout << endl; 
  print_m (lrsa, "lrsa");
  cout << endl << endl; 

  // part 2

  cm_t act (k, n); 
  rm_t art (k, n); 
  init_m (act, cls1()); 
  init_m (art, cls1());
  print_m (act, "act"); 
  cout << endl; 
  print_m (art, "art"); 
  cout << endl << endl;

  init_m (cmu, const_val<real_t> (0));
  init_m (cml, const_val<real_t> (0));
  init_m (rmu, const_val<real_t> (0));
  init_m (rml, const_val<real_t> (0));

  blas::syrk ( 1.0, bindings::trans(act), 0.0, bindings::upper(cmu)); 
  blas::syrk ( 1.0, bindings::trans(act), 0.0, bindings::lower(cml)); 
  blas::syrk ( 1.0, bindings::trans(art), 0.0, bindings::upper(rmu)); 
  blas::syrk ( 1.0, bindings::trans(art), 0.0, bindings::lower(rml)); 

  print_m (cmu, "cmu");
  cout << endl; 
  print_m (cml, "cml");
  cout << endl; 
  print_m (rmu, "rmu");
  cout << endl; 
  print_m (rml, "rml");
  cout << endl; 
}

