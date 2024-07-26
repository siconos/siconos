/* This file is adapted from (C) 2008 Maik Beckmann
 * https://lists.boost.org/MailArchives/ublas/2008/09/2984.php
 */

#ifndef determinant_hpp
#define determinant_hpp

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace siconos::externals::ublas {
template <class matrix_T>
double determinant(boost::numeric::ublas::matrix_expression<matrix_T> const& mat_r) {
  double det = 1.0;

  matrix_T mLu(mat_r());
  boost::numeric::ublas::permutation_matrix<std::size_t> pivots(mat_r().size1());

  int is_singular = lu_factorize(mLu, pivots);

  if (!is_singular) {
    for (std::size_t i = 0; i < pivots.size(); ++i) {
      if (pivots(i) != i) det *= -1.0;

      det *= mLu(i, i);
    }
  } else
    det = 0.0;

  return det;
}

}  // namespace siconos::externals::ublas
#endif
