//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "ublas_heev.hpp"

#include <boost/numeric/bindings/lapack/driver/hegv.hpp>
#include <boost/numeric/bindings/lapack/driver/sygv.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/hermitian.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/mpl/if.hpp>

#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;

template <typename T, typename W, typename UPLO>
int do_memory_uplo(int n, W& workspace)
{
  typedef typename bindings::remove_imaginary<T>::type real_type ;

  typedef ublas::matrix<T, ublas::column_major>     matrix_type ;
  typedef ublas::hermitian_adaptor<matrix_type, UPLO> hermitian_type ;
  typedef ublas::vector<real_type>                  vector_type ;

  // Set matrix
  matrix_type a(n, n);
  a.clear();
  vector_type e1(n);
  vector_type e2(n);

  fill(a);
  matrix_type a2(a);

  matrix_type b(n, n);
  b.clear();
  for(int i = 0; i < n; ++i) b(i,i) = 1;

  // Compute eigen decomposition.
  hermitian_type h_a(a);
  lapack::hegv(1, 'V', h_a, b, e1, workspace) ;

  if(check_residual(a2, e1, a)) return 255 ;

  hermitian_type h_a2(a2);
  lapack::hegv(1, 'N', h_a2, b, e2, workspace) ;
  if(norm_2(e1 - e2) > n * norm_2(e1) * std::numeric_limits< real_type >::epsilon()) return 255 ;

  // Test for a matrix range
  fill(a);
  a2.assign(a);

  typedef ublas::matrix_range< matrix_type > matrix_range ;
  typedef ublas::hermitian_adaptor<matrix_range, UPLO> hermitian_range_type ;

  ublas::range r(1,n-1) ;
  matrix_range a_r(a, r, r);
  ublas::vector_range< vector_type> e_r(e1, r);
  matrix_range b_r(b, r, r);

  hermitian_range_type h_a_r(a_r);
  lapack::hegv(1, 'V', h_a_r, b_r, e_r, workspace);

  matrix_range a2_r(a2, r, r);
  if(check_residual(a2_r, e_r, a_r)) return 255 ;

  return 0 ;
} // do_memory_uplo()


template <typename T, typename W>
int do_memory_type(int n, W workspace)
{
  std::cout << "  upper\n" ;
  if(do_memory_uplo<T,W,ublas::upper>(n, workspace)) return 255 ;
  std::cout << "  lower\n" ;
  if(do_memory_uplo<T,W,ublas::lower>(n, workspace)) return 255 ;
  return 0 ;
}


template <typename T>
struct Workspace
{
  typedef ublas::vector<T>                         array_type ;
  typedef lapack::detail::workspace1< array_type > type ;

  Workspace(size_t n)
    : work_(3*n-1)
  {}

  type operator()()
  {
    return lapack::workspace(work_) ;
  }

  array_type work_ ;
};

template <typename T>
struct Workspace< std::complex<T> >
{
  typedef ublas::vector<T>                                                 real_array_type ;
  typedef ublas::vector< std::complex<T> >                                 complex_array_type ;
  typedef lapack::detail::workspace2< complex_array_type,real_array_type > type ;

  Workspace(size_t n)
    : work_(2*n-1)
    , rwork_(3*n-2)
  {}

  type operator()()
  {
    return type(work_, rwork_);
  }

  complex_array_type work_ ;
  real_array_type    rwork_ ;
};

template <typename T>
int do_value_type()
{
  const int n = 8 ;

  std::cout << " optimal workspace\n";
  if(do_memory_type<T,lapack::optimal_workspace>(n, lapack::optimal_workspace())) return 255 ;

  std::cout << " minimal workspace\n";
  if(do_memory_type<T,lapack::minimal_workspace>(n, lapack::minimal_workspace())) return 255 ;

  std::cout << " workspace array\n";
  Workspace<T> work(n);
  if(do_memory_type<T,typename Workspace<T>::type >(n, work())) return 255 ;
  return 0;
} // do_value_type()


int main()
{
  // Run tests for different value_types
  std::cout << "float\n" ;
  if(do_value_type<float>()) return 255;

  std::cout << "double\n" ;
  if(do_value_type<double>()) return 255;

  std::cout << "complex<float>\n" ;
  if(do_value_type< std::complex<float> >()) return 255;

  std::cout << "complex<double>\n" ;
  if(do_value_type< std::complex<double> >()) return 255;

  std::cout << "Regression test succeeded\n" ;
  return 0;
}

