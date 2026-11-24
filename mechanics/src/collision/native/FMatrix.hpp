// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file FMatrix.hpp
 *  \brief Declaration of types used to handle plugged functions
 */

// Thif file must only contain forward declarations !  No implementation, no definitions !

#ifndef FMATRIX_DECL

#include <boost/numeric/ublas/fwd.hpp>  // boost::numeric fwd
#include <vector>

#include "FunctionTypes.hpp"

namespace siconos::collision::native {

/** A matrix of time function */
using FMatrix = boost::numeric::ublas::matrix<
    siconos::modeling::func_prototypes::FunctionS_S, boost::numeric::ublas::column_major,
    std::vector<siconos::modeling::func_prototypes::FunctionS_S>>;

}  // namespace siconos::collision::native

#endif
