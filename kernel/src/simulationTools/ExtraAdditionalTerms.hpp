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

/*! \file ExtraAdditionalTerms.hpp
 * \brief base class for struct of functions adding optional integration terms
 */

#ifndef ExtraAdditionalTerms_hpp
#define ExtraAdditionalTerms_hpp

#include "SiconosFwd.hpp"
#include "SimulationTypeDef.hpp"

struct ExtraAdditionalTerms {

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(ExtraAdditionalTerms);

public:
  /** initialize elements in the graph for the computations
   * \param DSG0 the graph of DynamicalSystems
   * \param model the current Model
   */
  virtual void init(DynamicalSystemsGraph& DSG0, const NonSmoothDynamicalSystem& nsds, const TimeDiscretisation & td) = 0;

  /** add smooth term to xfree (like the control input, the error correction for an observer)
   * \param DSG0 the graph of DynamicalSystems
   * \param dsgVD a DynamicalSystem in the DS graph
   * \param h the current timestep
   * \param xfree the free state to modify
   */
  virtual void addSmoothTerms(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double h, SiconosVector& xfree) = 0;

  /** add contribution to JacRhs for instance if \f$\dot{x} = f(x) + g(x)u\f$
   * \param DSG0 the graph of DynamicalSystems
   * \param dsgVD a DynamicalSystem in the DS graph
   * \param h the current timestep
   * \param jacRhs the jacobian to modify
   */
  virtual void addJacobianRhsContribution(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double h, SiconosMatrix& jacRhs) = 0;

  /** Desctructor */
  virtual ~ExtraAdditionalTerms() {};
};

#endif
