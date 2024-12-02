/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file FMatrix.hpp
 *  \brief Declaration of types used to handle plugged functions
 */

// Thif file must only contain forward declarations !  No implementation, no definitions !

#ifndef FMATRIX_DECL

#include <boost/numeric/ublas/fwd.hpp>  // boost::numeric fwd
#include <vector>

#include "PluginTypes.hpp"  // For FTime

namespace siconos::collision::native {

/** A matrix of time function */
using FMatrix =
    boost::numeric::ublas::matrix<siconos::plugins::FTime, boost::numeric::ublas::column_major,
                                  std::vector<siconos::plugins::FTime>>;

}  // namespace siconos::collision::native

#endif
