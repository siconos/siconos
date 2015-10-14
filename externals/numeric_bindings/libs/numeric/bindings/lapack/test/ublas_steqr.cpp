//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "../../blas/test/random.hpp"

#include <boost/numeric/bindings/lapack/computational/steqr.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <algorithm>
#include <limits>
#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;


template <typename T>
int do_value_type() {
   const int n = 10 ;

   typedef typename bindings::remove_imaginary<T>::type real_type ;
   typedef std::complex< real_type >                                            complex_type ;

   typedef ublas::matrix<T, ublas::column_major> matrix_type ;
   typedef ublas::vector<T>                      vector_type ;
   real_type safety_factor (1.5);

   // Set matrix
   matrix_type z( n, n );
   vector_type d( n ), e ( n - 1 ) ;

   std::fill( d.begin(), d.end(), 2.0 ) ;
   std::fill( e.begin(), e.end(), -1.0 ) ;

   // Compute eigendecomposition.
   lapack::steqr( 'I', d, e, z ) ;

   for ( int i=0; i<d.size(); ++i) {
     T sum( 0.0 ) ;
     for (int j=0; j<d.size(); ++j) {
       sum += z(i,j)*z(i,j) * d(j) ;
     }
     if (std::abs( sum - 2.0 ) > safety_factor*10 * std::numeric_limits<T>::epsilon() ) return 1 ;

     if (i>0) {
       sum = 0.0 ;
       for (int j=0; j<d.size(); ++j) {
         sum += z(i-1,j)*z(i,j) * d(j) ;
       }
       if (std::abs( sum + 1.0 ) > safety_factor*10 * std::numeric_limits<T>::epsilon() ) return 1 ;
     }
   }

   return 0 ;
} // do_value_type()



int main() {
   // Run tests for different value_types
   std::cout << "float\n" ;
   if (do_value_type<float>()) return 255;
   std::cout << "double\n" ;
   if (do_value_type<double>()) return 255;

   std::cout << "Regression test succeeded\n" ;
   return 0;
}

