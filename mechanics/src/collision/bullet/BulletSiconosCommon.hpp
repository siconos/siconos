/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#ifndef BulletSiconosCommon_hpp
#define BulletSiconosCommon_hpp

#include "SiconosVector.hpp"
#include <boost/math/quaternion.hpp>
#include <LinearMath/btVector3.h>
//void copyQuatRot(boost::math::quaternion<double>& from, SiconosVector& to);

void copyQuatRot(const SiconosVector& from, boost::math::quaternion<double>& to);

void copyQuatPos(const boost::math::quaternion<double>& from, SiconosVector& to);

void copyQuatPos(const SiconosVector& from, boost::math::quaternion<double>& to);

void copyQuatPos(const btVector3& from, boost::math::quaternion<double>& to);

void copyBtVector3(const btVector3 &from, SiconosVector& to);


//void copyQuatRot2d(boost::math::quaternion<double>& from, SiconosVector& to);


void copyQuatRot2d(const SiconosVector& from, boost::math::quaternion<double>& to);


void copyQuatPos2d(const boost::math::quaternion<double>& from, SiconosVector& to);


void copyQuatPos2d(const SiconosVector& from, boost::math::quaternion<double>& to);



void copyQuatPos2d(const btVector3& from, boost::math::quaternion<double>& to);

void copyBtVector32d(const btVector3 &from, SiconosVector& to);

void display_quat(boost::math::quaternion<double>& quat);
#endif
