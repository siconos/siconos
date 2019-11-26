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

#include "SimulationGraphs.hpp"
#include "ControlZOHAdditionalTerms.hpp"
#include "SimpleMatrixFriends.hpp"
#include "SiconosVectorFriends.hpp"

#include "Topology.hpp"
#include "MatrixIntegrator.hpp"
#include "SimpleMatrix.hpp"

//#define DEBUG_WHERE_MESSAGES

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

#include "debug.h"

void ControlZOHAdditionalTerms::init(DynamicalSystemsGraph& DSG0,
                                     const NonSmoothDynamicalSystem & nsds,
                                     const TimeDiscretisation & td)
{
  DEBUG_BEGIN("void ControlZOHAdditionalTerms::init(...)\n")
  DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  for (std11::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi)
  {
    DynamicalSystem& ds = *DSG0.bundle(*dsvi);
    if (DSG0.B.hasKey(*dsvi))
    {
      DSG0.Bd[*dsvi].reset(new MatrixIntegrator(ds, nsds, td, DSG0.B[*dsvi]));
      if (DSG0.Bd.at(*dsvi)->isConst())
        DSG0.Bd.at(*dsvi)->integrate();
    }
    if (DSG0.L.hasKey(*dsvi))
    {
      DSG0.Ld[*dsvi].reset(new MatrixIntegrator(ds, nsds, td, DSG0.L[*dsvi]));
      if (DSG0.Ld.at(*dsvi)->isConst())
        DSG0.Ld.at(*dsvi)->integrate();
    }
    if (DSG0.pluginB.hasKey(*dsvi))
      DSG0.Bd[*dsvi].reset(new MatrixIntegrator(ds, nsds, td, DSG0.pluginB[*dsvi], DSG0.u[*dsvi]->size()));
    if (DSG0.pluginL.hasKey(*dsvi))
      DSG0.Ld[*dsvi].reset(new MatrixIntegrator(ds, nsds, td, DSG0.pluginL[*dsvi], DSG0.e[*dsvi]->size()));
  }
  DEBUG_END("void ControlZOHAdditionalTerms::init(...)\n")
}

void ControlZOHAdditionalTerms::addSmoothTerms(DynamicalSystemsGraph& DSG0,
                                               const DynamicalSystemsGraph::VDescriptor& dsgVD,
                                               const double h, SiconosVector& xfree)
{
  DEBUG_BEGIN("void ControlZOHAdditionalTerms::addSmoothTerms(...)\n")
  // check whether we have a system with a control input
  if (DSG0.u.hasKey(dsgVD))
  {
    DEBUG_PRINT("a system has a control input\n");
    assert(DSG0.Bd.hasKey(dsgVD));
    if (!DSG0.Bd.at(dsgVD)->isConst())
    {
      DSG0.Bd.at(dsgVD)->integrate();
    }
    DEBUG_EXPR(DSG0.Bd.at(dsgVD)->mat().display());
    DEBUG_EXPR(xfree.display());
    DEBUG_EXPR((*DSG0.u.at(dsgVD)).display());
    prod(DSG0.Bd.at(dsgVD)->mat(), *DSG0.u.at(dsgVD), xfree, false); // xfree += Bd*u
  }
  // check whether the DynamicalSystem is an Observer
  if (DSG0.e.hasKey(dsgVD))
  {
    assert(DSG0.Ld.hasKey(dsgVD));
    if (!DSG0.Ld.at(dsgVD)->isConst())
    {
      DSG0.Ld.at(dsgVD)->integrate();
    }
    prod(DSG0.Ld.at(dsgVD)->mat(), *DSG0.e.at(dsgVD), xfree, false); // xfree += -Ld*e
  }
  DEBUG_END("void ControlZOHAdditionalTerms::addSmoothTerms(...)\n");
}

void ControlZOHAdditionalTerms::addJacobianRhsContribution(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double h, SiconosMatrix& jacRhs)
{
  // nothing to be done here
}
