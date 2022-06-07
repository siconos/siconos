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


// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

#include "BulletSiconosCommon.hpp"


#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif


#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic pop
#endif

#include <boost/math/quaternion.hpp>


//void copyQuatRot(boost::math::quaternion<double>& from, SiconosVector& to)
//{
//  to(3) = from.R_component_1();
//  to(4) = from.R_component_2();
//  to(5) = from.R_component_3();
//  to(6) = from.R_component_4();
//}

void copyQuatRot(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(from(3), from(4), from(5), from(6));
}


void copyQuatPos(const boost::math::quaternion<double>& from, SiconosVector& to)
{
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
  to(2) = from.R_component_4();
}


void copyQuatPos(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from(0), from(1), from(2));
}


void copyQuatPos(const btVector3& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from.x(), from.y(), from.z());
}

void copyBtVector3(const btVector3 &from, SiconosVector& to)
{
  to(0) = from.x();
  to(1) = from.y();
  to(2) = from.z();
}


//void copyQuatRot2d(boost::math::quaternion<double>& from, SiconosVector& to)
//{
//  double angle = 2.0 * acos(from.R_component_1());
//  to(2) = angle;
//}


void copyQuatRot2d(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  double half_angle = from(2)/2.0;

  to = boost::math::quaternion<double>(cos(half_angle), 0.0, 0.0, sin(half_angle));
}


void copyQuatPos2d(const boost::math::quaternion<double>& from, SiconosVector& to)
{
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
}


void copyQuatPos2d(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from(0), from(1), 0.0);
}


void copyQuatPos2d(const btVector3& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from.x(), from.y(), 0.0);
}

void copyBtVector32d(const btVector3 &from, SiconosVector& to)
{
  to(0) = from.x();
  to(1) = from.y();
}

void display_quat(boost::math::quaternion<double>& quat)
{
  std::cout << "q_0: " << quat.R_component_1()
            << " q_1: " << quat.R_component_2()
            << " q_2: " << quat.R_component_3()
            << " q_3: " << quat.R_component_4() << std::endl;
}
