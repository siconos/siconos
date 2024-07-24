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

#ifndef MBTB_CONTACTRELATION
#define MBTB_CONTACTRELATION

#include "MechanismsFwd.hpp"
#include <SiconosKernel.hpp>
#include "MBTB_Contact.hpp"
//! It is a relation dedicated for the simple unilateral (ie: without  Coulomb friction).
/*!
  Aggregation to the class MBTB_Contact, the member _pContact contains the CAD information.
  It derivates from siconos::NewtonEuler1DR. This class does the link between CAD and Siconos.
 */
class MBTB_ContactRelation : public NewtonEuler1DR
{

protected:
  MBTB_Contact * _pContact;
  MBTB_ContactRelation();
public:
  /** Constructor
   *  \param pC [in] a pointer to the MBTB_Contact. Must be allocated/free by the caller.
   */
  MBTB_ContactRelation(MBTB_Contact * pC);
  /** This function has to compute the distance between the objects.
   * \param time  the given time
   * \param q0 the position
   * \param y the output
   */
  virtual void computeh(double time, const BlockVector& q0, SiconosVector& y);
  //! Doing nothing.
  virtual ~MBTB_ContactRelation();

  ACCEPT_STD_VISITORS();

};


#endif
