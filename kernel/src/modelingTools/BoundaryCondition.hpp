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
#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "SiconosPointers.hpp"

#include "SiconosFwd.hpp"

#include "SiconosVector.hpp"
#include "PluggedObject.hpp"
#include "Tools.hpp"

typedef  void (*FPtrPrescribedVelocity)(double, unsigned int, double*);

/** \class BoundaryCondition
 *  \brief This class models simple boundary conditions for
 *   prescribing the velocities in a Dynamical System. A simple
 *   boundary condition is considered to fix a component \f$ j \f$ of
 *   * the velocity vector, i.e., \f$ v_j(t) = bc(t)\f$ where \f$
 *   bc(t)\f$ is a given function of time.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.2.0.
 *  \date (Re-Creation) September 2010
 *
 */
class BoundaryCondition
{
public:

  /** \fn BoundaryCondition(SP::UnsignedIntVector  newVelocityIndices);
   *  \brief Basic constructor
   *  \param newVelocityIndices the indices of the velocity subjected to prescribed velocities
   */

  BoundaryCondition(SP::UnsignedIntVector newVelocityIndices);

  /** \fn BoundaryCondition(SP::UnsignedIntVector  newVelocityIndices,
   *                 SP::SiconosVector newVelocityValues);
   *  \brief Constructor with constant prescribed values
   *  \param newVelocityIndices the indices of the velocity subjected to prescribed velocities
   *  \param newVelocityValues the values of the prescribed velocities
   */
  BoundaryCondition(SP::UnsignedIntVector  newVelocityIndices,
                    SP::SiconosVector newVelocityValues);

  /** destructor */
  virtual ~BoundaryCondition();

  // === GETTERS AND SETTERS ===

  /** to get the velocityIndices
   *  \return a pointer on _velocityIndices
   */
  inline SP::UnsignedIntVector velocityIndices()
  {
    return _velocityIndices;
  };

  /** to get the prescribedVelocity
   *  \return a pointer on _prescribedVelocity
   */
  inline SP::SiconosVector prescribedVelocity()
  {
    return _prescribedVelocity;
  };

  /** to get the prescribedVelocityOld
   *  \return a pointer on _prescribedVelocityOld
   */
  inline SP::SiconosVector prescribedVelocityOld()
  {
    return _prescribedVelocityOld;
  };

  /** allow to set a specified function to compute prescribedVelocity
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the name of the function to use in this plugin
   */
  void setComputePrescribedVelocityFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginPrescribedVelocity->setComputeFunction(pluginPath, functionName);
    if (!_prescribedVelocity) _prescribedVelocity.reset(new SiconosVector((unsigned int)_velocityIndices->size()));
  }

  /** default function to compute the precribed velocities
   *  \param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BoundaryCondition);

  /** protected default constructor */
  BoundaryCondition() {};

  /* Indices of the prescribed component of the velocity vector */
  SP::UnsignedIntVector _velocityIndices;

  /* Values of the prescribed component of the velocity vector */
  SP::SiconosVector _prescribedVelocity;

  /* Old values of the prescribed component of the velocity vector */
  SP::SiconosVector _prescribedVelocityOld;

  /*plugin defining the function V(t)*/
  SP::PluggedObject _pluginPrescribedVelocity;

  //   /*Link to the precribed DynamicalSystem*/
  //   SP::DynamicalSystem _DS;
};

TYPEDEF_SPTR(BoundaryCondition)
#endif // BOUNDARYCONDITION_HPP
