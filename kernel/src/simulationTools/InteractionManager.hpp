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

/*! \file InteractionManager.hpp
  \brief Definition of a class that manages dynamic interactions.
*/

#ifndef InteractionManager_h
#define InteractionManager_h

#include "Interaction.hpp"
#include "Model.hpp"
#include "DynamicalSystem.hpp"
#include "SiconosVisitor.hpp"
#include "NSLawMatrix.hpp"

class InteractionManager
{
public:
  InteractionManager() : _nslaws(1) {}
  virtual ~InteractionManager() {}

  /** Called by Simulation after updating positions prior to starting
   * the Newton loop. */
  virtual void updateInteractions(SP::Simulation simulation) {}

  /** Called by Simulation after updating positions inside the Newton
   * loop but prior to performing the next iteration.  Only override
   * if different behaviour during the loop is desired.  May or may
   * not be called depending on Simulation options. */
  virtual void updateInteractionsNewtonIteration(SP::Simulation simulation) {
    updateInteractions(simulation);
  }

  /** Specify a non-smooth law to use for a given combination of
   *  interaction groups.
   * \param nsl a SP::NonSmoothLaw
   * \param group1 first group
   * \param group2 second group */
  virtual void insertNonSmoothLaw(SP::NonSmoothLaw nslaw,
                                  long unsigned int group1,
                                  long unsigned int group2);

  /** Retrieve a non-smooth law to use for a given combination of
   *  interaction groups.
   * \return nsl a SP::NonSmoothLaw
   * \param group1 first group
   * \param group2 second group */
  virtual SP::NonSmoothLaw nonSmoothLaw(long unsigned int group1,
                                        long unsigned int group2);

protected:

  /** nslaws */
  NSLawMatrix _nslaws;

  friend class Simulation;
  friend class TimeStepping;
  friend class EventDriven;

  ACCEPT_SERIALIZATION(InteractionManager);
};

#endif /* InteractionManager_h */
