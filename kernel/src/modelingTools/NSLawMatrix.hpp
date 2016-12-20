/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/** \file NSLawMatrix.hpp
 *  \brief Base (abstract) class for a symmetric matrix of nonsmooth laws
*/

#ifndef NSLAWMATRIX_H
#define NSLAWMATRIX_H

#include <SiconosSerialization.hpp>
#include <NonSmoothLaw.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

/** This class uses (extends) a symmetric matrix to store
 *  NonSmoothLaws associated to pairs of integers.  It can be used
 *  e.g. to maintain a list of non-smooth laws associated with contact
 *  between types of objects. */
class NSLawMatrix : public boost::numeric::ublas::symmetric_matrix < SP::NonSmoothLaw >
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NSLawMatrix);

public:
  NSLawMatrix(NSLawMatrix::size_type i = 1)
    : boost::numeric::ublas::symmetric_matrix < SP::NonSmoothLaw >(i) {}
};

#endif
