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

/*! \file SphereLDSPlanR.hpp
  \brief SphereLDS relation with a plan - Inherits from LagrangianScleronomousR
*/

#ifndef SphereLDSPlanR_h
#define SphereLDSPlanR_h

#include "MechanicsFwd.hpp"
#include "CircularR.hpp"

class SphereLDSPlanR : public LagrangianScleronomousR, public std11::enable_shared_from_this<SphereLDSPlanR>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SphereLDSPlanR);


  /* Ax + By + Cz + D = 0 */
  double r, A, B, C, D, nN, nU;
  /* u ^ v  = n */
  double u1, u2, u3, v1, v2, v3, n1, n2, n3, ru1, ru2, ru3, rv1, rv2, rv3;

  SphereLDSPlanR() {};

public:

  /** Constructor
  \param r disk radius
  \param A
  \param B
  \param C
  \param D
  */
  SphereLDSPlanR(double r, double A, double B, double C, double D);

  double distance(double, double, double, double);

  using LagrangianScleronomousR::computeh;
  void computeh(SiconosVector& q, SiconosVector& z, SiconosVector& y);

  void computeJachq(SiconosVector& q, SiconosVector& z);

  bool equal(double _A, double _B, double _C, double _D, double _r) const
  {
    return (A == _A && B == _B && C == _C && D == _D && r == _r) ;
  }

  /** visitors hook
   */
  ACCEPT_VISITORS();

};
#endif /* SphereLDSPlanR_h */
