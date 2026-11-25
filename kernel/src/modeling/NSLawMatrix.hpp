/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/** \file NSLawMatrix.hpp
 *  \brief Base (abstract) class for a symmetric matrix of nonsmooth laws
 */

#ifndef NSLAWMATRIX_H
#define NSLAWMATRIX_H

#include <NonSmoothLaw.hpp>
#include <SiconosSerialization.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

namespace siconos::modeling {
/** This class uses (extends) a symmetric matrix to store
 *  NonSmoothLaws associated to pairs of integers.  It can be used
 *  e.g. to maintain a list of non-smooth laws associated with contact
 *  between types of objects. */
class NSLawMatrix
    : public boost::numeric::ublas::symmetric_matrix<std::shared_ptr<NonSmoothLaw>> {
 protected:
  ACCEPT_SERIALIZATION(NSLawMatrix);

 public:
  NSLawMatrix(NSLawMatrix::size_type i = 1)
      : boost::numeric::ublas::symmetric_matrix<std::shared_ptr<NonSmoothLaw>>(i) {}
};
}  // namespace siconos::modeling
#endif
