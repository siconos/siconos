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

/*! \file SphereLDS.hpp
  \brief Definition of a 3D sphere as a LagrangianDS (with Euler
         Angles)
*/


#ifndef SphereLDS_h
#define SphereLDS_h

#include "MechanicsFwd.hpp"
#include "LagrangianDS.hpp"

class SphereLDS : public LagrangianDS, public std11::enable_shared_from_this<SphereLDS>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SphereLDS);

  double radius;
  double massValue;
  double I;

  SphereLDS() {};

public:

  SphereLDS(double, double, SP::SiconosVector, SP::SiconosVector);

  ~SphereLDS();

  inline double getQ(unsigned int pos)
  {
    assert(pos < _ndof);
    return (*_q[0])(pos);
  };
  inline double getVelocity(unsigned int pos)
  {
    assert(pos < _ndof);
    return (*_q[1])(pos);
  };

  inline double getMassValue() const
  {
    return massValue;
  };

  inline double getRadius() const
  {
    return radius;
  };

  void computeMass();
  void computeMass(SP::SiconosVector)
  {
    RuntimeException::selfThrow("SphereLDS::computeMass(vector) - not implemented");
  }

  void computeFGyr(SP::SiconosVector, SP::SiconosVector);

  void computeFGyr();

  void computeJacobianFGyrq();
  void computeJacobianFGyrqDot();

  void computeJacobianFGyrq(SP::SiconosVector, SP::SiconosVector);
  void computeJacobianFGyrqDot(SP::SiconosVector, SP::SiconosVector);


  /** visitors hook
   */
  ACCEPT_BASE_VISITORS(LagrangianDS);

};
#endif /* SphereLDS_h */
