
// BLAS level 1 (vectors) 

//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <iostream>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/blas/level1.hpp>
#ifdef F_USE_STD_VECTOR
#include <vector>
#include <boost/numeric/bindings/std/vector.hpp> 
#endif 
#include "utils.h"

namespace blas = boost::numeric::bindings::blas;
namespace ublas = boost::numeric::ublas;

using std::cout;
using std::endl; 
using std::size_t; 

typedef double real_t;

#ifndef F_USE_STD_VECTOR
typedef ublas::vector<real_t> vct_t;
#else
typedef ublas::vector<real_t, std::vector<double> > vct_t;
#endif 

int main() {

  cout << endl; 

  vct_t v (10);
  init_v (v, times_plus<real_t> (0.1, 0.1)); 
  print_v (v, "v"); 
  vct_t vy (10); 
  blas::set (1., vy); 
  print_v (vy, "vy"); 
  cout << endl; 

  // v <- 2 v
  blas::scal (2.0, v); 
  print_v (v, "v <- 2 v");

  // vy <- 0.5 v + vy
  blas::axpy (0.5, v, vy); 
  print_v (vy, "vy <- 0.5 v + vy"); 

  // v^T vy
  cout << "v vy = " << blas::dot (v, vy) << endl;
  cout << endl; 

  /////////////////
  // ranges 

  // v -- new init 
  init_v (v, times_plus<real_t> (0.1, 0.1)); 
  print_v (v, "v"); 
  // vy -- new init 
  init_v (vy, kpp (1)); 
  print_v (vy, "vy"); 

  // v[2..6]
  ublas::vector_range<vct_t> vr (v, ublas::range (2, 6)); 
  print_v (vr, "v[2..6]"); 
  // vy[4..8] 
  ublas::vector_range<vct_t> vry (vy, ublas::range (4, 8)); 
  print_v (vry, "vy[4..8]"); 

  // v[2..6] <- 0.1 v[2..6]
  blas::scal (0.1, vr); 
  print_v (v, "v[2..6] <- 0.1 v[2..6]"); 

  // vr^T vr 
  // ublas::vector_range<vct_t const> cvr (v, ublas::range (2, 6)); 
  cout << "v[2..6] v[2..6] = " << blas::dot (vr, vr) << endl;

  // vy[4..8] <- v[2..6] + vy[4..8]
  blas::axpy (1.0, vr, vry); 
  print_v (vy, "vy[4..8] <- v[2..6] + vy[4..8]"); 
  cout << endl; 

  /////////////////
  // slices 

  // v -- new init 
  init_v (v, times_plus<real_t> (0.1, 0.1)); 
  print_v (v, "v"); 
  // vy -- new init 
  init_v (vy, kpp (1)); 
  print_v (vy, "vy"); 

  // v[1:2:4]
  ublas::vector_slice<vct_t> vs (v, ublas::slice (1, 2, 4)); 
  print_v (vs, "v[1:2:4]"); 
  // vy[2:2:4] 
  ublas::vector_slice<vct_t> vsy (vy, ublas::slice (2, 2, 4)); 
  print_v (vsy, "vy[2:2:4]"); 

  // v[1:2:4] <- 10 v[1:2:4]
  blas::scal (10.0, vs); 
  print_v (v, "v[1:2:4] <- 10 v[1:2:4]"); 

  // vs^T vs
  cout << "v[1:2:4] v[1:2:4] = " << blas::dot (vs, vs) << endl;

  // vy[2:2:4] <- 0.01 v[1:2:4] + vy[2:2:4] 
  blas::axpy (0.01, vs, vsy); 
  print_v (vy, "vy[2:2:4] <- 0.01 v[1:2:4] + vy[2:2:4]"); 
  cout << endl; 

  ////////////////////////////////////////////
  // ranges & slices 

  // v -- new init 
  init_v (v, times_plus<real_t> (0.1, 0.1)); 
  print_v (v, "v"); 
  // vy <- 1.0
  blas::set (1., vy);
  print_v (vy, "vy <- 1.0"); 

  // vy[2:2:4] <- 0.01 v[2..6] + vy[2:2:4] 
  blas::axpy (0.01, vr, vsy); 
  print_v (vy, "vy[2:2:4] <- 0.01 v[2..6] + vy[2:2:4]"); 

  // vy <- 1.0
  blas::set (1., vy); 
  print_v (vy, "vy <- 1.0"); 

  // vy[4..8] <- 0.01 v[1:2:4] + vy[4..8] 
  blas::axpy (0.01, vs, vry); 
  print_v (vy, "vy[4..8] <- 0.01 v[1:2:4] + vy[4..8]"); 

  // vr^T vs == vs^T vr
  cout << "v[2..6] v[1:2:4] = " << blas::dot (vr, vs) << " == " 
       << blas::dot (vs, vr) << endl;

  // vr^T vsy == vsy^T vr
  cout << "v[2..6] vy[2:2:4] = " << blas::dot (vr, vsy) << " == " 
       << blas::dot (vsy, vr) << endl;

  cout << endl; 

  ///////////////////
  // slice of range

  // v -- new init 
  init_v (v, times_plus<real_t> (0.1, 0.1)); 
  print_v (v, "v"); 

  // v[1..9][1:2:3] 
  ublas::vector_range<vct_t> vr1 (v, ublas::range (1, 9)); 
  ublas::vector_slice< 
    ublas::vector_range<vct_t> 
  > vsr (vr1, ublas::slice (1, 2, 3)); 
  print_v (vr1, "v[1..9]"); 
  print_v (vsr, "v[1..9][1:2:3]"); 

  // v[1..9][1:2:3] <- 0.1 v[1..9][1:2:3]
  blas::scal (0.1, vsr); 
  print_v (v, "0.1 v[1..9][1:2:3]"); 

  // ||vsr||_2 
  cout << "||v[1..9][1:2:3]||_2 = " << blas::nrm2 (vsr) << endl; 
  cout << endl; 

  ///////////////////
  // range of slice 

  // vy <- 1.0
  init_v (vy, kpp (1)); 
  print_v (vy, "vy"); 

  // v[0:2:5][1:4] 
  ublas::vector_slice<vct_t> vsy1 (vy, ublas::slice (0, 2, 5));
  ublas::vector_range<
    ublas::vector_slice<vct_t> 
  > vrs (vsy1, ublas::range (1, 4)); 
  print_v (vsy1, "v[0:2:5]");
  print_v (vrs, "v[0:2:5][1:4]");

  // v[0:2:5][1:4] <- 0.01 v[1..9][1:2:3] + v[0:2:5][1:4]
  blas::axpy (0.01, vsr, vrs); 
  print_v (vy, "0.01 v[1..9][1:2:3] + v[0:2:5][1:4]"); 

  cout << endl; 

}
