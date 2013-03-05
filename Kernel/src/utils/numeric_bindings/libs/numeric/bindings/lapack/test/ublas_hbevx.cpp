//
//  Copyright Toon Knapen, Karl Meerbergen
//  Copyright Thomas Klimpel 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "ublas_heev.hpp"

#include <boost/numeric/bindings/lapack/driver/hbevx.hpp>
#include <boost/numeric/bindings/lapack/driver/sbevx.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/banded.hpp>
#include <boost/numeric/bindings/ublas/hermitian.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/mpl/if.hpp>

#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;

template <typename U>
int lower() {
   return 0;
}

template <>
int lower<ublas::lower>() { return 2; }


template <typename U>
int upper() {
   return 0;
}

template <>
int upper<ublas::upper>() { return 2; }


template <typename T, typename W, typename UPLO, typename Orientation>
int do_memory_uplo(int n, W& workspace ) {
   typedef typename bindings::remove_imaginary<T>::type real_type ;

   typedef ublas::banded_matrix<T, Orientation>      banded_type ;
   typedef ublas::matrix<T, ublas::column_major>     matrix_type ;
   typedef ublas::vector<real_type>                  vector_type ;

   // Set matrix
   banded_type a( n, n, lower<UPLO>(), upper<UPLO>() ); a.clear(); // Make upper triangular band matrix
   matrix_type q( n, n );
   matrix_type z( n, n );
   vector_type e1( n );
   vector_type e2( n );

   fill_banded( a );
   banded_type a2( a );

   ublas::hermitian_adaptor<banded_type, UPLO> h( a ), h2( a2 );
   real_type vl = 0, vu = 1, abstol = 1e-28;
   fortran_int_t il = 1, iu = n-1, m;
   ublas::vector<fortran_int_t> ifail(n);


   // Compute Schur decomposition.
   lapack::hbevx( 'V', 'A',
     h, q, vl, vu, il, iu, abstol, m, e1, z, ifail, workspace ) ;

   if (check_residual( h2, e1, z )) return 255 ;

   lapack::hbevx( 'N', 'A',
     h2, q, vl, vu, il, iu, abstol, m, e2, z, ifail, workspace ) ;
   if (norm_2( e1 - e2 ) > n * norm_2( e1 ) * std::numeric_limits< real_type >::epsilon()) return 255 ;

   // Test for a matrix range
   fill_banded( a ); a2.assign( a );

   typedef ublas::matrix_range< banded_type > banded_range ;

   ublas::range r(1,n-1) ;
   banded_range a_r( a, r, r );
   ublas::hermitian_adaptor< banded_range, UPLO> h_r( a_r );
   ublas::vector_range< vector_type> e_r( e1, r );
   ublas::matrix_range< matrix_type> z_r( z, r, r );
   ublas::vector<fortran_int_t> ifail_r(n-2);

   lapack::hbevx( 'V', 'A',
     h_r, q, vl, vu, il, iu, abstol, m, e_r, z_r, ifail_r, workspace ) ;

   banded_range a2_r( a2, r, r );
   ublas::hermitian_adaptor< banded_range, UPLO> h2_r( a2_r );
   if (check_residual( h2_r, e_r, z_r )) return 255 ;

   return 0 ;
} // do_memory_uplo()


template <typename T, typename W>
int do_memory_type(int n, W workspace) {
   std::cout << "  upper\n" ;
   if (do_memory_uplo<T,W,ublas::upper, ublas::row_major>(n, workspace)) return 255 ;
   std::cout << "  lower\n" ;
   if (do_memory_uplo<T,W,ublas::lower, ublas::row_major>(n, workspace)) return 255 ;
   return 0 ;
}

template <typename T>
struct Workspace {
   typedef ublas::vector< T >                                      array_type ;
   typedef ublas::vector< fortran_int_t >                              int_array_type ;
   typedef lapack::detail::workspace2<array_type, int_array_type > type ;

   Workspace(size_t n)
   : work_( 7*n )
   , iwork_( 5*n )
   {}

   type operator() () {
      return lapack::workspace(work_ , iwork_) ;
   }

   array_type work_ ;
   int_array_type iwork_ ;
};

template <typename T>
struct Workspace< std::complex<T> > {
   typedef ublas::vector< std::complex<T> >                 complex_array_type ;
   typedef ublas::vector< T >                               real_array_type ;
   typedef ublas::vector< fortran_int_t >                       int_array_type ;
   typedef lapack::detail::workspace3<
      complex_array_type, real_array_type, int_array_type > type ;

   Workspace(size_t n)
   : work_( n )
   , rwork_( 7*n )
   , iwork_( 5*n )
   {}

   type operator() () {
      return lapack::workspace(work_, rwork_, iwork_) ;
   }

   complex_array_type work_ ;
   real_array_type    rwork_ ;
   int_array_type     iwork_ ;
};

template <typename T>
int do_value_type() {
   const int n = 8 ;

   std::cout << " optimal workspace\n";
   if (do_memory_type<T,lapack::optimal_workspace>( n, lapack::optimal_workspace() ) ) return 255 ;

   std::cout << " minimal workspace\n";
   if (do_memory_type<T,lapack::minimal_workspace>( n, lapack::minimal_workspace() ) ) return 255 ;

   std::cout << " workspace array\n";
   Workspace<T> work( n );
   if (do_memory_type<T,typename Workspace<T>::type >( n, work() ) ) return 255 ;
   return 0;
} // do_value_type()


int main() {
   // Run tests for different value_types
   std::cout << "float\n" ;
   if (do_value_type<float>()) return 255;

   std::cout << "double\n" ;
   if (do_value_type<double>()) return 255;

   std::cout << "complex<float>\n" ;
   if (do_value_type< std::complex<float> >()) return 255;

   std::cout << "complex<double>\n" ;
   if (do_value_type< std::complex<double> >()) return 255;

   std::cout << "Regression test succeeded\n" ;
   return 0;
}

