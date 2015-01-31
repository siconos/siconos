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
#include "ControlZOHAdditionalTerms.hpp"

#include "Topology.hpp"
#include "MatrixIntegrator.hpp"

void ControlZOHAdditionalTerms::init(DynamicalSystemsGraph& DSG0, const Model& model)
{
  DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  for (std11::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi)
  {
    DynamicalSystem& ds = *DSG0.bundle(*dsvi);
    if (DSG0.B.hasKey(*dsvi))
    {
      DSG0.Bd[*dsvi].reset(new MatrixIntegrator(ds, model, DSG0.B[*dsvi]));
      if (DSG0.Bd.at(*dsvi)->isConst())
        DSG0.Bd.at(*dsvi)->integrate();
    }
    if (DSG0.L.hasKey(*dsvi))
    {
      DSG0.Ld[*dsvi].reset(new MatrixIntegrator(ds, model, DSG0.L[*dsvi]));
      if (DSG0.Ld.at(*dsvi)->isConst())
        DSG0.Ld.at(*dsvi)->integrate();
    }
    if (DSG0.pluginB.hasKey(*dsvi))
      DSG0.Bd[*dsvi].reset(new MatrixIntegrator(ds, model, DSG0.pluginB[*dsvi], DSG0.u[*dsvi]->size()));
    if (DSG0.pluginL.hasKey(*dsvi))
      DSG0.Ld[*dsvi].reset(new MatrixIntegrator(ds, model, DSG0.pluginL[*dsvi], DSG0.e[*dsvi]->size()));
  }
}

void ControlZOHAdditionalTerms::addSmoothTerms(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double h, SiconosVector& xfree)
{
  // check whether we have a system with a control input
  if (DSG0.u.hasKey(dsgVD))
  {
    assert(DSG0.Bd.hasKey(dsgVD));
    if (!DSG0.Bd.at(dsgVD)->isConst())
    {
      DSG0.Bd.at(dsgVD)->integrate();
    }
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
}

void ControlZOHAdditionalTerms::addJacobianRhsContribution(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double h, SiconosMatrix& jacRhs)
{
  // nothing to be done here
}
