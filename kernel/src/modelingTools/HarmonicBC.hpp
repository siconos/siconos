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
#ifndef HARMONICBC_HPP
#define HARMONICBC_HPP


#include "BoundaryCondition.hpp"

/**\brief This class models a simple harmonic boundary conditions for
 *   prescribing the velocities in a Dynamical System. A simple
 *   boundary condition is considered to fix a component \f$ j \f$ of
 *   the velocity vector, i.e., \f$ v_j(t) = a +  b cos( \omega t+ \phi) \f$.
 *
 */
class HarmonicBC : public  BoundaryCondition
{
public:

  /** Constructor
   * \param newVelocityIndices the indices of the velocity subjected to prescribed velocities
   * \param a constant value for additive term of the prescribed velocity
   * \param b constant value for multiplicative term of the prescribed velocity
   * \param omega frequency
   * \param phi phase
   */
  HarmonicBC(SP::UnsignedIntVector newVelocityIndices,
             double a, double b,
             double omega, double phi) ;

  HarmonicBC(SP::UnsignedIntVector newVelocityIndices,
             SP::SiconosVector a, SP::SiconosVector b,
             SP::SiconosVector omega, SP::SiconosVector phi);


  /** destructor */
  virtual ~HarmonicBC();

  /** default function to compute the precribed velocities
   *
   *  \param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

protected:
  
  ACCEPT_SERIALIZATION(HarmonicBC);

  /** protected default constructor */
  HarmonicBC(): BoundaryCondition() {};

  /** Constant additive term of the prescribed velocity  */
  double _a;
  /** Constant multiplicative term of the prescribed velocity  */
  double _b;
  /** Constant frequency  */
  double _omega;
  /** Constant phase  */
  double _phi;

  /** Constant additive term of the prescribed velocity  */
  SP::SiconosVector _aV;
  /** Constant multiplicative term of the prescribed velocity  */
  SP::SiconosVector _bV;
  /** Constant frequency  */
  SP::SiconosVector _omegaV;
  /** Constant phase  */
  SP::SiconosVector _phiV;

  

  //   /*Link to the precribed DynamicalSystem*/
  //   SP::DynamicalSystem _DS;
};

TYPEDEF_SPTR(HarmonicBC)
#endif // HARMONICBC_HPP
