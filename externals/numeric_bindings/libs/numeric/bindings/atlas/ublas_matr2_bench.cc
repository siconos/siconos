
// BLAS level 2
// benchmarks 

//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

//#define USE_STD_VECTOR
//#define BOUNDED 100*100

//#define PRINT
//#define PRINT_M

//#define MODIFY

#include <stddef.h>
#include <iostream>
#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include <boost/numeric/ublas/io.hpp>
#ifdef USE_STD_VECTOR
#include <vector>
#include <boost/numeric/bindings/std/vector.hpp> 
#endif 
#include <boost/timer.hpp>
#include "utils.h" 

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;

using std::cout;
using std::cin;
using std::endl; 

typedef double real_t; 

#ifdef USE_STD_VECTOR
typedef std::vector<real_t> storage_t; 
#else
#ifndef BOUNDED 
typedef ublas::unbounded_array<real_t> storage_t; 
#else
typedef ublas::bounded_array<real_t, BOUNDED> storage_t; 
#endif 
#endif 

typedef ublas::vector<real_t, storage_t> vct_t;

typedef ublas::matrix<real_t, ublas::column_major, storage_t> cm_t;
typedef ublas::matrix<real_t, ublas::row_major, storage_t> rm_t;

typedef ublas::symmetric_adaptor<cm_t, ublas::upper> ucsa_t; 
typedef ublas::symmetric_adaptor<cm_t, ublas::lower> lcsa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::upper> ursa_t; 
typedef ublas::symmetric_adaptor<rm_t, ublas::lower> lrsa_t; 

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

////////////////////////////////////////////////////
// general matrix: gemv()

template <typename M>
void bench_gemv (size_t n, size_t runs, const char* msg) {

  cout << msg << endl; 

  vct_t x (n);
  blas::set (1., x);
  vct_t y (n); 

  M a (n, n);
  init_symm (a); 
#ifdef PRINT_M
  print_m (a, "a");
  cout << endl; 
#endif 

  boost::timer (t); 
  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    a (0, 2) = r; 
    a (2, 0) = r; 
#endif 
    blas::gemv ( 1.0, a, x, 0.0, y);
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif
  } 
  cout << "  gemv:        " << t.elapsed() << endl;  

  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    a (0, 2) = r; 
    a (2, 0) = r; 
#endif 
    blas::gemv ( 1., bindings::trans(a), x, 0., y );
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif
  } 
  cout << "  gemv trans:  " << t.elapsed() << endl;  

  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    a (0, 2) = r; 
    a (2, 0) = r; 
#endif 
    y = prod (a, x);
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif 
  }
  cout << "  = prod:      " << t.elapsed() << endl; 

  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    a (0, 2) = r; 
    a (2, 0) = r; 
#endif 
    y.assign (prod (a, x));
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif
  } 
  cout << "  assign prod: " << t.elapsed() << endl; 
  cout << endl; 
}


/////////////////////////////////////////////////////////
// symmetric adaptor: symv()

template <typename M, typename SA>
void bench_symv (size_t n, size_t runs, const char* msg) {

  cout << msg << endl; 

  vct_t x (n);
  blas::set (1., x);
  vct_t y (n); 

  M a (n, n);
  SA sa (a); 
  init_symm (sa, 'l'); 
#ifdef PRINT_M
  print_m (sa, "sa"); 
  cout << endl; 
  print_m (a, "a"); 
  cout << endl; 
#endif 

  boost::timer (t); 
  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    sa (0, 2) = r; 
#endif 
    blas::symv ( 1., sa, x, 0., y);
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif
  } 
  cout << "  symv:        " << t.elapsed() << endl;  

  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    sa (0, 2) = r; 
#endif 
    y = prod (sa, x);
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif 
  }
  cout << "  = prod:      " << t.elapsed() << endl; 

  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    sa (0, 2) = r; 
#endif 
    y.assign (prod (sa, x));
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif
  } 
  cout << "  assign prod: " << t.elapsed() << endl; 
  cout << endl; 
}


/////////////////////////////////////////////////////////
// symmetric matrix: spmv()

template <typename SM>
void bench_spmv (size_t n, size_t runs, const char* msg) {

  cout << msg << endl; 

  vct_t x (n);
  blas::set (1., x);
  vct_t y (n); 

  SM sa (n, n);
  init_symm (sa, 'l'); 
#ifdef PRINT_M
  cout << sa << endl;
  cout << endl; 
#endif 

  boost::timer (t); 
  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    sa (0, 2) = r; 
#endif 
    blas::spmv ( 1.0, sa, x, 0.0, y);
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif
  } 
  cout << "  spmv:        " << t.elapsed() << endl;  

  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    sa (0, 2) = r; 
#endif 
    y = prod (sa, x);
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif 
  }
  cout << "  = prod:      " << t.elapsed() << endl; 

  t.restart(); 
  for (size_t r = 0; r < runs; ++r) {
#ifdef MODIFY
    sa (0, 2) = r; 
#endif 
    y.assign (prod (sa, x));
#ifdef PRINT
    std::cout << "y " << bindings::noop( y ) << std::endl;
#endif
  } 
  cout << "  assign prod: " << t.elapsed() << endl; 
  cout << endl; 
}


//////////////////////////////////////////////////////
int main (int argc, char **argv) {
  size_t n = 0, r = 0;
  if (argc > 1) {
    n = atoi(argv [1]);
  }
  if (argc > 2) {
    r = atoi(argv [2]);
  }

  cout << endl; 

  if (n <= 0) {
    cout << "n -> ";
    cin >> n;
  }
  if (r <= 0) {
    cout << "r -> ";
    cin >> r;
  }
  cout << endl; 

  bench_gemv<rm_t> (n, r, "row major"); 
  bench_gemv<cm_t> (n, r, "column major"); 

  bench_symv<rm_t, ursa_t> (n, r, "symmetric ad, row, upper"); 
  bench_symv<rm_t, lrsa_t> (n, r, "symmetric ad, row, lower"); 
  bench_symv<cm_t, ucsa_t> (n, r, "symmetric ad, column, upper"); 
  bench_symv<cm_t, lcsa_t> (n, r, "symmetric ad, column, lower"); 

  bench_spmv<ursymm_t> (n, r, "symmetric, row, upper"); 
  bench_spmv<lrsymm_t> (n, r, "symmetric, row, lower"); 
  bench_spmv<ucsymm_t> (n, r, "symmetric, column, upper"); 
  bench_spmv<lcsymm_t> (n, r, "symmetric, column, lower"); 

}
