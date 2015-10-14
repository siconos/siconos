
// BLAS level 2 -- complex numbers

#include <iostream>
#include <complex>
#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/conj.hpp>
#include <boost/numeric/bindings/io.hpp>
#include <boost/numeric/bindings/noop.hpp>
#ifdef F_USE_STD_VECTOR
#include <boost/numeric/bindings/std/vector.hpp> 
#endif 
#include "utils.h" 

namespace ublas = boost::numeric::ublas;
namespace blas = boost::numeric::bindings::blas;

using std::cout;
using std::endl; 

typedef double real_t;
typedef std::complex<real_t> cmplx_t; 

#ifndef F_USE_STD_VECTOR
typedef ublas::vector<cmplx_t> vct_t;
typedef ublas::matrix<cmplx_t, ublas::row_major> m_t;
#else
typedef ublas::vector<cmplx_t, std::vector<cmplx_t> > vct_t;
typedef ublas::matrix<cmplx_t, ublas::column_major, std::vector<cmplx_t> > m_t;
#endif 

int main() {

  cout << endl; 

  vct_t vx (2);
  blas::set( 1, vx );
  std::cout << "vx " << bindings::noop( vx ) << std::endl;
  vct_t vy (4); // vector size can be larger 
                // than corresponding matrix size 
  blas::set( 0, vy );
  std::cout << "vy " << bindings::noop( vy ) << std::endl;
  cout << endl; 

  m_t m (3, 2);
  init_m (m, kpp (1)); 
  print_m (m, "m"); 
  cout << endl; 

  // vy = m vx
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  cout << endl; 

  m (0, 0) = cmplx_t (0., 1.);
  m (0, 1) = cmplx_t (0., 2.);
  m (1, 0) = cmplx_t (0., 3.);
  m (1, 1) = cmplx_t (0., 4.);
  m (2, 0) = cmplx_t (0., 5.);
  m (2, 1) = cmplx_t (0., 6.);
  print_m (m, "m"); 
  cout << endl; 

  // vy = m vx
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  cout << endl; 

  m (0, 0) = cmplx_t (-1., 1.);
  m (0, 1) = cmplx_t (-2., 2.);
  m (1, 0) = cmplx_t (-3., 3.);
  m (1, 1) = cmplx_t (-4., 4.);
  m (2, 0) = cmplx_t (-5., 5.);
  m (2, 1) = cmplx_t (-6., 6.);
  print_m (m, "m"); 
  cout << endl; 

  // vy = m vx
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  blas::gemv ( 1.0, m, vx, 0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  cout << endl; 

  blas::set ( 1, vx );
  std::cout << "vx " << bindings::noop( vx ) << std::endl;

  // vy = m vx
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  cout << endl; 

  blas::set ( cmplx_t (1, 1), vx );
  std::cout << "vx " << bindings::noop( vx ) << std::endl;

  // vy = m vx
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  blas::gemv ( 1.0, m, vx, 0.0, vy);
  std::cout << "m vx " << bindings::noop( vy ) << std::endl;
  cout << endl; 

  // vx = m^H vy
  blas::set ( cmplx_t (-1,-1), vy );
  std::cout << "vy " << bindings::noop( vy ) << std::endl;
  blas::gemv ( 1.0, bindings::conj(m), vy, 0.0, vx);
  std::cout << "m^H vy " << bindings::noop( vx ) << std::endl;
  cout << endl; 

  m_t mx (2, 2); 
  m_t my (3, 2); 


  ublas::matrix_column<m_t> mxc0 (mx, 0), 
                            mxc1 (mx, 1); 
  ublas::matrix_column<m_t> myc0 (my, 0),
                            myc1 (my, 1); 

  blas::set ( cmplx_t (1, 0), mxc0 );
  blas::set ( cmplx_t (0, 0), mxc1 );
  blas::set ( cmplx_t (0, 0), myc0 );
  blas::set ( cmplx_t (0, 0), myc1 );

  print_m (mx, "mx");
  cout << endl; 
  print_m (my, "my");
  cout << endl; 

  // my[.,0] = m mx[.,0] 
  blas::gemv ( 1.0, m, mxc0, 0.0, myc0); 
  print_m (my, "m mx[.,0]");

  cout << endl;

}
