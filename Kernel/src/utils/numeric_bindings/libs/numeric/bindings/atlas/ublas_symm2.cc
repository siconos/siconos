
// BLAS level 2
// symmetric & hermitian matrices

//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <stddef.h>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/blas/level1/set.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/ublas/hermitian.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;
namespace bindings = boost::numeric::bindings;

using std::cout;
using std::cin;
using std::endl; 

typedef double real_t; 
typedef std::complex<real_t> cmplx_t; 

typedef ublas::vector<real_t> vct_t;
typedef ublas::vector<cmplx_t> cvct_t;

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

typedef ublas::hermitian_matrix<
  cmplx_t, ublas::upper, ublas::column_major
> ucherm_t; 
typedef ublas::hermitian_matrix<
  cmplx_t, ublas::lower, ublas::column_major
> lcherm_t; 
typedef ublas::hermitian_matrix<
  cmplx_t, ublas::upper, ublas::row_major
> urherm_t; 
typedef ublas::hermitian_matrix<
  cmplx_t, ublas::lower, ublas::row_major
> lrherm_t; 

typedef ublas::matrix<cmplx_t, ublas::column_major> cm_t; 
typedef ublas::matrix<cmplx_t, ublas::row_major> rm_t; 

typedef ublas::hermitian_adaptor<cm_t, ublas::upper> ucha_t; 
typedef ublas::hermitian_adaptor<cm_t, ublas::lower> lcha_t; 
typedef ublas::hermitian_adaptor<rm_t, ublas::upper> urha_t; 
typedef ublas::hermitian_adaptor<rm_t, ublas::lower> lrha_t; 


int main (int argc, char **argv) {
  size_t n = 0;
  if (argc > 1) {
    n = atoi(argv [1]);
  }
  
  cout << endl; 

  cout << "symmetric matrix" << endl << endl; 
  if (n <= 0) {
    cout << "n -> ";
    cin >> n;
  }
  cout << endl; 

  vct_t vx (n), vy (n); 
  blas::set( 1., vx );
  print_v (vx, "vx"); 

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
  print_m_data (ucs, "ucs");
  cout << endl; 

  print_m (lcs, "lcs");
  cout << endl; 
  print_m_data (lcs, "lcs");
  cout << endl; 

  print_m (urs, "urs");
  cout << endl; 
  print_m_data (urs, "urs");
  cout << endl; 

  print_m (lrs, "lrs");
  cout << endl; 
  print_m_data (lrs, "lrs");
  cout << endl; 

  // vy = symm vx 
  blas::spmv( 1, ucs, vx, 0, vy ); 
  print_v (vy, "vy = ucs vx"); 
  cout << endl; 
  blas::spmv (1, lcs, vx, 0, vy); 
  print_v (vy, "vy = lcs vx"); 
  cout << endl; 
  blas::spmv (1., urs, vx, 0., vy); 
  print_v (vy, "vy = urs vx"); 
  cout << endl; 
  blas::spmv (1.0f, lrs, vx, 0.0f, vy); 
  print_v (vy, "vy = lrs vx"); 
  cout << endl; 

#ifdef F_COMPILATION_FAILURE

  blas::symv (ucs, vx, vy); 
  print_v (vy, "vy = ucs vx"); 
  cout << endl; 

  blas::hpmv (lcs, vx, vy); 
  print_v (vy, "vy = lcs vx"); 
  cout << endl; 

  blas::gemv (urs, vx, vy); 
  print_v (vy, "vy = urs vx"); 
  cout << endl; 

  blas::spmv (cmplx_t (1., 0.), lrs, vx, cmplx_t (0., 0.), vy); 
  print_v (vy, "vy = lrs vx"); 
  cout << endl; 

#endif 

  ///////////////////////////////////////////////////

  cout << endl << "hermitian matrix" << endl << endl; 

  size_t n2 = 3; 

  cvct_t cvx (n2), cvy (n2); 
  blas::set( 1., cvx );
  print_v (cvx, "cvx"); 
  cout << endl; 

  ucherm_t uch (n2, n2); 
  uch (0, 0) = cmplx_t (3, 0); 
  uch (0, 1) = cmplx_t (2, -2); 
  uch (1, 1) = cmplx_t (3, 0); 
  uch (0, 2) = cmplx_t (1, -1); 
  uch (1, 2) = cmplx_t (2, 2); 
  uch (2, 2) = cmplx_t (3, 0); 
  print_m (uch, "uch"); 
  cout << endl; 
  print_m_data (uch, "uch"); 
  cout << endl; 

  lcherm_t lch (n2, n2); 
  lch = uch; 
  print_m (lch, "lch"); 
  cout << endl; 
  print_m_data (lch, "lch"); 
  cout << endl; 

  urherm_t urh (uch); 
  print_m (urh, "urh"); 
  cout << endl; 
  print_m_data (urh, "urh"); 
  cout << endl; 

  lrherm_t lrh (uch); 
  print_m (lrh, "lrh"); 
  cout << endl; 
  print_m_data (lrh, "lrh"); 
  cout << endl; 

  // cvy = herm cvx 
  blas::hpmv (1, uch, cvx, 0, cvy); 
  print_v (cvy, "cvy = uch cvx"); 
  cout << endl; 
  blas::hpmv (1, lch, cvx, 0, cvy); 
  print_v (cvy, "cvy = lch cvx"); 
  cout << endl; 
  blas::hpmv (1., urh, cvx, 0., cvy); 
  print_v (cvy, "cvy = urh cvx"); 
  cout << endl; 
  blas::hpmv (cmplx_t (1., 0.), lrh, cvx, cmplx_t (0., 0.), cvy); 
  print_v (cvy, "cvy = lrh cvx"); 
  cout << endl; 

  ///////////////////////////////////////////////////

  cout << endl << "hermitian adaptor" << endl << endl; 

  cm_t cm1 (n2, n2); 
  ucha_t ucha (cm1);
  ucha = uch; 
  print_m (ucha, "ucha"); 
  cout << endl; 
  print_m_data (ucha, "ucha"); 
  cout << endl; 

  cm_t cm2 (n2, n2); 
  lcha_t lcha (cm2);
  lcha = uch;  
  print_m (lcha, "lcha"); 
  cout << endl; 
  print_m_data (lcha, "lcha"); 
  cout << endl; 

  rm_t rm1 (n2, n2); 
  urha_t urha (rm1); 
  urha = uch; 
  print_m (urha, "urha"); 
  cout << endl; 
  print_m_data (urha, "urha"); 
  cout << endl; 

  rm_t rm2 (n2, n2); 
  lrha_t lrha (rm2); 
  lrha = uch; 
  print_m (lrha, "lrha"); 
  cout << endl; 
  print_m_data (lrha, "lrha"); 
  cout << endl; 

  // cvy = herma cvx 
  blas::hemv (1, ucha, cvx, 0, cvy); 
  print_v (cvy, "cvy = uch cvx"); 
  cout << endl; 
  blas::hemv (1, lcha, cvx, 0, cvy); 
  print_v (cvy, "cvy = lch cvx"); 
  cout << endl; 
  blas::hemv (1., urha, cvx, 0., cvy); 
  print_v (cvy, "cvy = urh cvx"); 
  cout << endl; 
  blas::hemv (cmplx_t (1., 0.), lrha, cvx, cmplx_t (0., 0.), cvy); 
  print_v (cvy, "cvy = lrh cvx"); 
  cout << endl; 

}
