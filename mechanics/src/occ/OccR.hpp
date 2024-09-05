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
/** \file OccR.hpp
    \brief A Siconos Newton Euler 3D friction relation between
    two BRep contact points.
 */

#ifndef OccR_hpp
#define OccR_hpp

#include "Geometer.hpp"        // For  DistanceCalculator
#include "NewtonEuler3DR.hpp"  // Base class

namespace siconos::mechanics::occ {

struct ContactPoint;

class OccR : public siconos::modeling::NewtonEuler3DR {
 protected:
  const ContactPoint& _contact1;
  const ContactPoint& _contact2;

  DistanceCalculator dt{CadmbtbDistanceType{}};

  //  std::shared_ptr<Geometer> _geometer{nullptr};

  double _offset1{0.};
  double _offset2{0.};

 public:
  /** Constructor from contact points.
   *
   *  \param contact1 : the first contact.
   *  \param contact2 : the second contact.
   */
  OccR(const ContactPoint& contact1, const ContactPoint& contact2,
       const DistanceCalculator& distance_calculator);

  // Rule of 5
  OccR() = delete;
  OccR(const OccR&) = delete;
  OccR(OccR&&) = delete;
  OccR& operator=(const OccR&) = delete;
  OccR& operator=(OccR&&) = delete;
  ~OccR() noexcept = default;

  /** Compute h.
   *
   *  \param time : the time.
   *  \param q0 : the state vector.
   *  \param y : output vector.
   */
  void computeh(double time, const siconos::algebra::BlockVector& q0,
                siconos::algebra::SiconosVector& y);

  /** Set offset1, offset from first contact.
   *
   *  \param val : the new value.
   */
  void setOffset1(double val) { _offset1 = val; };

  /** Set offset2, offset from second contact.
   *
   *  \param val : the new value.
   */
  void setOffset2(double val) { _offset2 = val; };

  // /** Get geometer.
  //  *
  //  *  \return a std::shared_ptr<Geometer> object.
  //  */
  // std::shared_ptr<Geometer> geometer() { return _geometer; }

  // /** Set geometer.
  //  *
  //  *  \param geometer the new geometer
  //  */
  // void setGeometer(std::shared_ptr<Geometer> geometer) { _geometer = geometer; }
};
}  // namespace siconos::mechanics::occ
#endif
