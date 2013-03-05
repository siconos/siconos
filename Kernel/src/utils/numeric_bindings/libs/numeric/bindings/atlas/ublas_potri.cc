
// inverting symmetric/hermitian positive definite
// factor (potrf()) and invert (potri())

// #define BOOST_UBLAS_STRICT_HERMITIAN
// .. doesn't work (yet?)  

//#define BOOST_NO_FUNCTION_TEMPLATE_ORDERING

#include <cstddef>
#include <iostream>
#include <complex>
#include <boost/numeric/bindings/blas/level3.hpp>
#include <boost/numeric/bindings/lapack/computational/potri.hpp>
#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/ublas/hermitian.hpp>
#include "utils.h"

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;
namespace bindings = boost::numeric::bindings;
namespace lapack = boost::numeric::bindings::lapack;

using std::size_t; 
using std::cout;
using std::endl; 

typedef double real_t; 

typedef std::complex<real_t> cmplx_t; 

#ifndef F_ROW_MAJOR
typedef ublas::matrix<real_t, ublas::column_major> m_t;
typedef ublas::matrix<cmplx_t, ublas::column_major> cm_t;
#else
typedef ublas::matrix<real_t, ublas::row_major> m_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> cm_t;
#endif

#ifndef F_UPPER
typedef ublas::symmetric_adaptor<m_t, ublas::lower> symm_t; 
typedef ublas::hermitian_adaptor<cm_t, ublas::lower> herm_t; 
#else
typedef ublas::symmetric_adaptor<m_t, ublas::upper> symm_t; 
typedef ublas::hermitian_adaptor<cm_t, ublas::upper> herm_t; 
#endif 

int main() {

  cout << endl; 

  cout << "real symmetric\n" << endl; 

  size_t n = 3; 
  m_t a (n, n);    // matrix (storage)
  symm_t sa (a);   // symmetric adaptor 

#ifdef F_UPPER
  init_symm (sa, 'u'); 
#else
  init_symm (sa, 'l'); 
#endif
  // ifdef F_UPPER 
  //        [5 4 3 2 1]
  //        [0 5 4 3 2]
  //    a = [0 0 5 4 3]
  //        [0 0 0 5 4]
  //        [0 0 0 0 n]
  // else 
  //        [5 0 0 0 0]
  //        [4 5 0 0 0]
  //    a = [3 4 5 0 0]
  //        [2 3 4 5 0]
  //        [1 2 3 4 5]
  print_m (sa, "A"); 
  cout << endl; 

  m_t a2 (sa);   // full symmetric copy of sa:
                 // .. sa is `lost' after potrf(); 
                 // .. only one parameter of symm() is symmetric matrix

  int ierr = lapack::potrf (sa); 
  if (!ierr) {
    lapack::potri (sa); 
    // ri should be (almost) identity matrix: 
    m_t ri (n, n); 
    blas::symm ( bindings::tag::right(), 1.0, sa, a2, 0.0, ri); 
    print_m (ri, "I = A * A^(-1)"); 
    cout << endl; 
    blas::symm ( bindings::tag::left(), 1.0, sa, a2, 0.0, ri); 
    print_m (ri, "I = A^(-1) * A"); 
    cout << endl; 
  }

  cout << "\n===========================\n" << endl; 
  cout << "complex hermitian (almost ;o)\n" << endl; 

  // hermitian 
  cm_t ca (3, 3); 
  herm_t ha (ca); 

#ifndef F_UPPER
  ha (0, 0) = cmplx_t (3, 0);
  ha (1, 0) = cmplx_t (2, 0);
  ha (1, 1) = cmplx_t (3, 0);
  ha (2, 0) = cmplx_t (1, 0);
  ha (2, 1) = cmplx_t (2, 0);
  ha (2, 2) = cmplx_t (3, 0);
#else
  ha (0, 0) = cmplx_t (3, 0);
  ha (0, 1) = cmplx_t (2, 0);
  ha (0, 2) = cmplx_t (1, 0);
  ha (1, 1) = cmplx_t (3, 0);
  ha (1, 2) = cmplx_t (2, 0);
  ha (2, 2) = cmplx_t (3, 0);
#endif 

  print_m (ha, "A"); 
  cout << endl; 

  cm_t ca2 (ha);  // full hermitian 
  
  ierr = lapack::potri (ha);   // potrf()
  if (ierr == 0) {
    lapack::potri (ha);        // potri()
    cm_t ic (3, 3); 
    blas::hemm ( bindings::tag::right(), 1.0, ha, ca2, 0.0, ic); 
    print_m (ic, "I = A * A^(-1)"); 
    cout << endl; 
  }
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl; 


  cout << "\n===========================\n" << endl; 
  cout << "complex hermitian\n" << endl; 

#ifndef F_UPPER
  ha (0, 0) = cmplx_t (25, 0);
  ha (1, 0) = cmplx_t (-5, 5);
  ha (1, 1) = cmplx_t (51, 0);
  ha (2, 0) = cmplx_t (10, -5);
  ha (2, 1) = cmplx_t (4, 6);
  ha (2, 2) = cmplx_t (71, 0);
#else
  ha (0, 0) = cmplx_t (25, 0);
  ha (0, 1) = cmplx_t (-5, -5);
  ha (0, 2) = cmplx_t (10, 5);
  ha (1, 1) = cmplx_t (51, 0);
  ha (1, 2) = cmplx_t (4, -6);
  ha (2, 2) = cmplx_t (71, 0);
#endif
  print_m (ha, "A"); 
  cout << endl; 

  ca2 = ha; 
  
  ierr = lapack::potrf (ha); 
  if (ierr == 0) {
    lapack::potri (ha); 
    cm_t ic (3, 3); 
    blas::hemm ( bindings::tag::right(), 1.0, ha, ca2, 0.0, ic); 
    print_m (ic, "I = A * A^(-1)"); 
    cout << endl; 
    blas::hemm ( bindings::tag::left(), 1.0, ha, ca2, 0.0, ic); 
    print_m (ic, "I = A^(-1) * A"); 
    cout << endl; 
  }
  else 
    cout << "matrix is not positive definite: ierr = " 
         << ierr << endl; 

  cout << endl; 

}

