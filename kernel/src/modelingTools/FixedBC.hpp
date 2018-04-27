/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#ifndef FIXEDBC_HPP
#define FIXEDBC_HPP


#include "BoundaryCondition.hpp"

/** \class FixedBC
 *  \brief This class models a simple fixed boundary conditions for
 *   prescribing the velocities in a Dynamical System. A simple
 *   boundary condition is considered to fix a component \f$ j \f$ of
 *   the velocity vector, i.e., \f$ v_j(t) = 0\f$ 
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 4.1.0
 *  \date November 2016
 *
 */
class FixedBC : public  BoundaryCondition
{
public:

  /** \fn FixedBC(SP::UnsignedIntVector  newVelocityIndices);
   *  \brief Basic constructor
   *  \param newVelocityIndices the indices of the velocity subjected to prescribed velocities
   */

  FixedBC(SP::UnsignedIntVector newVelocityIndices) ;

  /** destructor */
  virtual ~FixedBC();

  /** default function to compute the precribed velocities
   *  \param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FixedBC);

  /** protected default constructor */
  FixedBC(): BoundaryCondition() {};

};

TYPEDEF_SPTR(FixedBC)
#endif // FIXEDBC_HPP
