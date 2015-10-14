//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//


#include "../../blas/test/random.hpp"

#include <boost/numeric/bindings/lapack/computational/sytrd.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <algorithm>
#include <limits>
#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;


template <typename T, typename UPLO>
int do_value_type() {
   const int n = 10 ;

   typedef typename bindings::remove_imaginary<T>::type real_type ;
   typedef std::complex< real_type >                                            complex_type ;

   typedef ublas::matrix<T, ublas::column_major> matrix_type ;
   typedef ublas::symmetric_adaptor<matrix_type, UPLO> symmetric_type ;
   typedef ublas::vector<T>                      vector_type ;

   // Set matrix
   matrix_type a( n, n );
   vector_type d( n ), e ( n - 1 ), tau( n - 1 ) ;

   for (int i=0; i<n; ++i ) {
     for (int j=0; j<n; ++j ) {
       a(j,i) = 0.0 ;
     }
   }

   a(0,0) = 2.0 ;
   for (int i=1; i<n; ++i ) {
      a(i,i) = 2.0 ;
      a(i-1,i) = -1.0 ;
   }

   // Compute eigendecomposition.
   symmetric_type s_a( a );
   lapack::sytrd( s_a, d, e, tau ) ;

   for ( int i=0; i<d.size(); ++i) {
      if (std::abs( d(i) - 2.0 ) > 10 * std::numeric_limits<T>::epsilon() ) return 1 ;
   }
   for ( int i=0; i<e.size(); ++i) {
      if (std::abs( e(i) + 1.0 ) > 10 * std::numeric_limits<T>::epsilon() ) return 1 ;
   }

   return 0 ;
} // do_value_type()



int main() {
   // Run tests for different value_types
   if (do_value_type<float, ublas::upper>()) return 255;
   if (do_value_type<double, ublas::upper>()) return 255;

   std::cout << "Regression test succeeded\n" ;
   return 0;
}

