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

#ifndef MBTB_MOREAU_H
#define MBTB_MOREAU_H

#include "MoreauJeanOSI.hpp"  // Base class

namespace siconos::mechanisms {
/**
 * \brief This class implements a variant of the std MoreauJeanOSI TS
 * It inherits from Siconos::MoreauJeanOSI
 * the main variants lies in the activation and desactivation of constraints
 */
class MBTB_MoreauJeanOSI : public siconos::integrators::MoreauJeanOSI {
 public:
  double _deactivateYPosThreshold{1e-4};
  double _deactivateYVelThreshold{0.};
  double _activateYPosThreshold{0.};
  double _activateYVelThreshold{100.};

 public:
  /** constructor from a minimum set of data: one DS and its theta
   *  \param theta value for the theta parameter (default = 0.5)
   *  \param gamma value for the gamma parameter (default = NaN and gamma is not
   * used)
   */
  MBTB_MoreauJeanOSI(double theta = 0.5,
                     double gamma = std::numeric_limits<double>::quiet_NaN())
      : siconos::integrators::MoreauJeanOSI{theta, gamma} {};

  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   * \param inter interaction
   * \param i level
   * \return a Boolean
   */
  bool addInteractionInIndexSet(
      std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   * \param inter interaction
   * \param i level
   * \return a Boolean
   */
  bool removeInteractionFromIndexSet(
      std::shared_ptr<siconos::modeling::Interaction> inter, unsigned int i);
};
}  // namespace siconos::mechanisms
#endif
