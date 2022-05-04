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

/*! \file MultipleImpactNSL.hpp 
 */

#ifndef _MULTIPLEIMPACTNSL_
#define _MULTIPLEIMPACTNSL_

#include "NonSmoothLaw.hpp"
#include "SiconosPointers.hpp"      // for DEFINE_SPTR
#include "SiconosSerialization.hpp" // for ACCEPT_SERIALIZATION
#include "SiconosVisitor.hpp"       // for ACCEPT_STD_VISITORS


class MultipleImpactNSL : public NonSmoothLaw {
private:

  ACCEPT_SERIALIZATION(MultipleImpactNSL);

  // Energytical restitution coefficient
  double _ResCof;
  // Normal stiffness at contact
  double _Stiff;
  // Elasticity coefficient
  double _ElasCof;

public:
  // Default Constructor
  MultipleImpactNSL();
  // Constructor with parameters
  MultipleImpactNSL(double, double, double, unsigned int _dim = 1);
  // Destructor
  ~MultipleImpactNSL();
  // Get the value of the energytical restitution coefficientx
  inline double ResCof() const { return _ResCof; };
  // Get the value of the stiffness
  inline double Stiff() const { return _Stiff; };
  // Get the value of the elasticity coefficient
  inline double ElasCof() const { return _ElasCof; }
  // Set the value to the restitution coefficient
  void setResCof(double newResCof);
  // Set the value to the stiffness
  void setStiff(double newStiff);
  // Set the value to the elasticity cofficient
  void setElasCoeff(double _newElasCoef);
  //
  bool isVerified() const override;
  // Display the information about the multiple impact law
  void display() const override;
  // Visitors hook
  ACCEPT_STD_VISITORS();
};
DEFINE_SPTR(MultipleImpactNSL)
#endif
