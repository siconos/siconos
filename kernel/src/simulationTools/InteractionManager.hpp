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

/*! \file InteractionManager.hpp
  \brief Definition of a class that manages dynamic interactions.
*/

#ifndef InteractionManager_h
#define InteractionManager_h

#include <memory>

#include "NSLawMatrix.hpp"
#include "SiconosSerialization.hpp"

namespace siconos::modeling {
class NonSmoothLaw;
class NSLawMatrix;
}  // namespace siconos::modeling

namespace siconos::simulation {

class Simulation;

  class InteractionManager : public std::enable_shared_from_this<InteractionManager>{
  InteractionManager(const InteractionManager&) = delete;
  InteractionManager(InteractionManager&&) = delete;
  InteractionManager operator=(const InteractionManager&&) = delete;
  InteractionManager operator=(InteractionManager&&) = delete;

 public:
  InteractionManager() = default;

  virtual ~InteractionManager() noexcept = default;

  /** Called by Simulation after updating positions prior to starting
   * the Newton loop. */
  virtual void updateInteractions(std::shared_ptr<siconos::simulation::Simulation> simulation)
  {
  }

  /** Specify a non-smooth law to use for a given combination of
   *  interaction groups.
   * \param nslaw the new nonsmooth law
   * \param group1 id of the fisrt group of interactions
   * \param group2  id of the second group of interactions
   */
  virtual void insertNonSmoothLaw(std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw,
                                  unsigned long int group1, unsigned long int group2);

  /** Retrieve a non-smooth law to use for a given combination of
   *  interaction groups.
   * \return nsl a std::shared_ptr<siconos::modeling::NonSmoothLaw>
   * \param group1 first group
   * \param group2 second group */
  virtual std::shared_ptr<siconos::modeling::NonSmoothLaw> nonSmoothLaw(
      unsigned long int group1, unsigned long int group2);

 protected:
  /** nslaws */
  siconos::modeling::NSLawMatrix _nslaws{1};

  friend class Simulation;
  friend class TimeStepping;
  friend class EventDriven;

  ACCEPT_SERIALIZATION(InteractionManager);
};
}  // namespace siconos::simulation
#endif /* InteractionManager_h */
