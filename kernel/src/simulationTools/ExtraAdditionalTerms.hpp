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

/*! \file ExtraAdditionalTerms.hpp
 * \brief base class for struct of functions adding optional integration terms
 */

#ifndef ExtraAdditionalTerms_hpp
#define ExtraAdditionalTerms_hpp

#include "SimulationGraphs.hpp"

namespace siconos::modeling {
class NonSmoothDynamicalSystem;
}

namespace siconos::simulation {
class TimeDiscretisation;
}

namespace siconos::integrators {

/** Pure virtual class useful to define extra terms in the one-step integrators.*/
struct ExtraAdditionalTerms {
 private:
  ACCEPT_SERIALIZATION(ExtraAdditionalTerms);

 public:
  /** initialize elements in the graph for the computations
   *
   *  \param DSG0 the graph of DynamicalSystems
   *  \param nsds the current nonsmooth dynamical system
   *  \param td the current time discretisation
   */
  virtual void init(siconos::graphs::DynamicalSystemsGraph& DSG0,
                    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                    std::shared_ptr<siconos::simulation::TimeDiscretisation> td) = 0;

  /** add smooth term to xfree (like the control input, the error correction for an observer)
   *
   *  \param DSG0 the graph of DynamicalSystems
   *  \param dsgVD a DynamicalSystem in the DS graph
   *  \param h the current timestep
   *  \param xfree the free state to modify
   */
  virtual void addSmoothTerms(siconos::graphs::DynamicalSystemsGraph& DSG0,
                              const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD,
                              const double h, siconos::algebra::SiconosVector& xfree) = 0;

  /** add contribution to JacRhs for instance if \f$\dot{x} = f(x) + g(x)u\f$
   *
   *  \param DSG0 the graph of DynamicalSystems
   *  \param dsgVD a DynamicalSystem in the DS graph
   *  \param h the current timestep
   *  \param jacRhs the jacobian to modify
   */
  virtual void addJacobianRhsContribution(
      siconos::graphs::DynamicalSystemsGraph& DSG0,
      const siconos::graphs::DynamicalSystemsGraph::VDescriptor& dsgVD, const double h,
      siconos::algebra::SiconosMatrix& jacRhs) = 0;

  /** Desctructor */
  virtual ~ExtraAdditionalTerms() noexcept = default;
};
}  // namespace siconos::integrators
#endif
