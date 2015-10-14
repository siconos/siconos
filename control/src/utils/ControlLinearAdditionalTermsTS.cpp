/* Siconos-Kernel, Copyright INRIA 2005-2014
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include "SimulationGraphs.hpp"
#include "ControlLinearAdditionalTermsTS.hpp"

#include "Topology.hpp"
#include "MatrixIntegrator.hpp"

void ControlLinearAdditionalTermsTS::init(DynamicalSystemsGraph& DSG0, const Model& model)
{
  // Do nothing here
}


void ControlLinearAdditionalTermsTS::addSmoothTerms(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double h, SiconosVector& xfree)
{
  // check whether we have a system with a control input
  if (DSG0.u.hasKey(dsgVD))
  {
    assert(DSG0.B.hasKey(dsgVD));
    prod(h, *DSG0.B[dsgVD], *DSG0.u[dsgVD], xfree, false); // xfree += h*B*u
  }
  // check whether the DynamicalSystem is an Observer
  if (DSG0.e.hasKey(dsgVD))
  {
    assert(DSG0.L.hasKey(dsgVD));
    prod(h, *DSG0.L[dsgVD], *DSG0.e[dsgVD], xfree, false); // xfree += -h*L*e
  }
}

void ControlLinearAdditionalTermsTS::addJacobianRhsContribution(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double t, SiconosMatrix& jacRhs)
{
  // nothing to be done here
}
