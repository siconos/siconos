//
//  Copyright Toon Knapen, Karl Meerbergen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "../../blas/test/random.hpp"

#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/lapack/driver/gees.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/vector_view.hpp>
//#include <boost/numeric/bindings/detail/complex_utils.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/mpl/if.hpp>

#include <iostream>
#include <limits>


namespace ublas = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;
namespace bindings = boost::numeric::bindings;

struct apply_real {
  template< typename MatrixA, typename VectorW, typename MatrixVS,
        typename Workspace >
  static inline std::ptrdiff_t gees( const char jobvs, const char sort,
        external_fp select, MatrixA& a, fortran_int_t& sdim, VectorW& w,
        MatrixVS& vs, Workspace work ) {
    return lapack::gees( jobvs, sort, select, a, sdim, w, vs, work );
    /*
    fortran_int_t info = lapack::gees( jobvs, sort, select, a, sdim,
      bindings::detail::real_part_view(w), bindings::detail::imag_part_view(w),
      vs, work );
    bindings::detail::interlace(w);
    return info;
    */
  }
};

struct apply_complex {
  template< typename MatrixA, typename VectorW, typename MatrixVS,
        typename Workspace >
  static inline std::ptrdiff_t gees( const char jobvs, const char sort,
        external_fp select, MatrixA& a, fortran_int_t& sdim, VectorW& w,
        MatrixVS& vs, Workspace work ) {
    return lapack::gees( jobvs, sort, select, a, sdim, w, vs, work );
  }
};


// Randomize a matrix
template <typename M>
void randomize(M& m) {
   typedef typename M::size_type  size_type ;
   typedef typename M::value_type value_type ;

   size_type size1 = m.size1() ;
   size_type size2 = m.size2() ;

   for (size_type i=0; i<size2; ++i) {
      for (size_type j=0; j<size1; ++j) {
         m(j,i) = random_value< value_type >() ;
      }
   }
} // randomize()


template< typename T >
struct dispatch_select
{
};
template<>
struct dispatch_select<float>
{
  static fortran_bool_t my_select(float* w_real, float* w_imag) {
    return *w_real > std::abs(*w_imag);
  }
};
template<>
struct dispatch_select<double>
{
  static fortran_bool_t my_select(double* w_real, double* w_imag) {
    return *w_real > std::abs(*w_imag);
  }
};
template<>
struct dispatch_select<std::complex<float> >
{
  static fortran_bool_t my_select(std::complex<float>* w) {
    return w->real() > std::abs(w->imag());
  }
};
template<>
struct dispatch_select<std::complex<double> >
{
  static fortran_bool_t my_select(std::complex<double>* w) {
    return w->real() > std::abs(w->imag());
  }
};


template <typename T, typename W>
int do_memory_type(int n, W workspace) {
   typedef typename boost::mpl::if_<boost::is_complex<T>, apply_complex, apply_real>::type apply_t;
   typedef typename bindings::remove_imaginary< T>::type real_type ;
   typedef std::complex< real_type >                                            complex_type ;

   typedef ublas::matrix<T, ublas::column_major> matrix_type ;
   typedef ublas::vector<complex_type>           vector_type ;
   double safety_factor (1.5);

   // Set matrix
   matrix_type a( n, n );
   matrix_type z( n, n );
   vector_type e1( n );
   vector_type e2( n );

   randomize( a );
   matrix_type a2( a );
   external_fp select = reinterpret_cast<external_fp>(dispatch_select<T>::my_select);
   fortran_int_t sdim_info(0);
   // Compute Schur decomposition.
   apply_t::gees( 'V', 'S', select, a, sdim_info, e1, z, workspace ) ;

   // Check Schur factorization
   if (norm_frobenius( prod( a2, z ) - prod( z, a ) )
           >= safety_factor*10.0* norm_frobenius( a2 ) * std::numeric_limits< real_type >::epsilon() ) return 255 ;

   matrix_type z_dummy( 1, 1 );
   apply_t::gees( 'N', 'S', select, a2, sdim_info, e2, z_dummy, workspace ) ;
   if (norm_2( e1 - e2 ) > safety_factor*norm_2( e1 ) * std::numeric_limits< real_type >::epsilon()) return 255 ;

   if (norm_frobenius( a2 - a )
           >= safety_factor*10.0* norm_frobenius( a2 ) * std::numeric_limits< real_type >::epsilon() ) return 255 ;


   return 0 ;
} // do_value_type()


template <typename T>
struct Workspace {
   typedef ublas::vector<T>                         array_type ;
   typedef ublas::vector< fortran_bool_t >               bool_array_type ;
   typedef lapack::detail::workspace2< array_type,bool_array_type > type ;

   Workspace(size_t n)
   : work_( 3*n )
   , bwork_( n )
   {}

   type operator() () {
      return lapack::workspace(work_, bwork_) ;
   }

   array_type work_ ;
   bool_array_type bwork_;
};


template <typename T>
struct Workspace< std::complex<T> > {
   typedef ublas::vector<T>                                                 real_array_type ;
   typedef ublas::vector< std::complex<T> >                                 complex_array_type ;
   typedef ublas::vector< fortran_bool_t >                                       bool_array_type ;
   typedef lapack::detail::workspace3< complex_array_type,real_array_type,bool_array_type > type ;

   Workspace(size_t n)
   : work_( 2*n )
   , rwork_( n )
   , bwork_( n )
   {}

   type operator() () {
      return lapack::workspace(work_, rwork_, bwork_) ;
   }

   complex_array_type work_ ;
   real_array_type    rwork_ ;
   bool_array_type    bwork_;
};


template <typename T>
int do_value_type() {
   const int n = 8 ;
   
   if (do_memory_type<T,lapack::optimal_workspace>( n, lapack::optimal_workspace() ) ) return 255 ;
   if (do_memory_type<T,lapack::minimal_workspace>( n, lapack::minimal_workspace() ) ) return 255 ;

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
