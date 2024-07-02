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

#ifndef MBTB_FC3DCONTACTRELATION
#define MBTB_FC3DCONTACTRELATION

#include "NewtonEuler3DR.hpp"  // Base class

namespace siconos::mechanisms {

class MBTB_Contact;

//! It is a relation dedicated for the unilateral constraint with Coulomb
//! friction.
/*!
  Aggregation to the class MBTB_Contact, the member _pContact contains the CAD
  information. It derivates from Siconos::NewtonEuler3DR. This class does the
  link between CAD and Siconos.
 */

class MBTB_FC3DContactRelation : public siconos::modeling::NewtonEuler3DR {
 protected:
  //! The member containing the CAD properties.
  std::shared_ptr<MBTB_Contact> _pContact{nullptr};

  // MBTB_FC3DContactRelation();

 public:
  /** Constructor
   * \param _pContact [in] a pointer to the MBTB_Contact. Must be allocated/free
   * by the caller.
   */
  MBTB_FC3DContactRelation(std::shared_ptr<MBTB_Contact> pContact);

  /** This function has to compute the distance between the objects.
   * \param time the given  time
   * \param q0 the position
   * \param y the output
   */
  virtual void computeh(double time, const siconos::algebra::BlockVector& q0,
                        siconos::algebra::SiconosVector& y) override;

  virtual ~MBTB_FC3DContactRelation() noexcept = default;
};
}  // namespace siconos::mechanisms
#endif
