/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/** \file BulletSiconosCommon.hpp
    Library of function to convert between siconos or boost quaternions and bullet vector
    Note: these functions are in addition of those of RotationQuaternion.hpp (kernel)
*/

#ifndef BulletSiconosCommon_hpp
#define BulletSiconosCommon_hpp

#include <boost/math_fwd.hpp>  // for quaternion

#include "BulletDeclarations.h"  // for btVector3
#include "SiconosVector.hpp"

namespace siconos::collision::bullet {

void copyQuatPos(const btVector3& from, boost::math::quaternion<double>& to);

void copyBtVector3(const btVector3& from, siconos::algebra::SiconosVector3& to);

void copyQuatPos2d(const btVector3& from, boost::math::quaternion<double>& to);

void copyBtVector32d(const btVector3& from, siconos::algebra::SiconosVector2& to);

}  // namespace siconos::collision::bullet
#endif
