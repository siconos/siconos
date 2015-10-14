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
#include "ControlLinearAdditionalTermsED.hpp"

#include "Topology.hpp"
#include "MatrixIntegrator.hpp"

typedef void (*AdditionalTermsEDfctU)(double, unsigned, double*, unsigned, double*, double*, unsigned, double*);

void ControlLinearAdditionalTermsED::init(DynamicalSystemsGraph& DSG0, const Model& model)
{

  DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  for (std11::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi)
  {
    DynamicalSystem& ds = *DSG0.bundle(*dsvi);
    if (DSG0.pluginU.hasKey(*dsvi))
    {
      DSG0.tmpXdot[*dsvi].reset(new SiconosVector(ds.getx().size()));
    }
    if (DSG0.pluginJacgx.hasKey(*dsvi))
    {
      DSG0.jacgx[*dsvi].reset(new SimpleMatrix(ds.getx().size(), ds.getx().size()));
    }

  }
}


void ControlLinearAdditionalTermsED::addSmoothTerms(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double t, SiconosVector& xdot)
{
  // check whether we have a system with a control input
  if (DSG0.u.hasKey(dsgVD))
  {
    if (DSG0.B.hasKey(dsgVD))
    {
      prod(DSG0.B.getRef(dsgVD), DSG0.u.getRef(dsgVD), xdot, false); // xdot += B*u
    }
    else if (DSG0.pluginU.hasKey(dsgVD))
    {
      DynamicalSystem& ds = *DSG0.bundle(dsgVD);
      SiconosVector& u = DSG0.u.getRef(dsgVD);
      SiconosVector& tmpXdot = DSG0.tmpXdot.getRef(dsgVD);
      ((AdditionalTermsEDfctU)DSG0.pluginU.getRef(dsgVD).fPtr)(t, xdot.size(), ds.getx().getArray(), u.size(), u.getArray(), tmpXdot.getArray(), ds.getz().size(), ds.getz().getArray());
      xdot += tmpXdot; // xdot += g(x, u)
    }
    else
    {
      RuntimeException::selfThrow("ControlLinearAdditionalTermsED :: input u but no B nor pluginU");
    }
  }
  // check whether the DynamicalSystem is an Observer
  if (DSG0.e.hasKey(dsgVD))
  {
    assert(DSG0.L.hasKey(dsgVD));
    prod(*DSG0.L[dsgVD], *DSG0.e[dsgVD], xdot, false); // xdot += -L*e
  }
}

void ControlLinearAdditionalTermsED::addJacobianRhsContribution(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double t, SiconosMatrix& jacRhs)
{
  // check whether we have a system with a control input
  if (DSG0.pluginJacgx.hasKey(dsgVD))
  {
    DynamicalSystem& ds = *DSG0.bundle(dsgVD);
    SiconosVector& u = DSG0.u.getRef(dsgVD);
    SimpleMatrix& tmpJacgx = DSG0.jacgx.getRef(dsgVD);
    ((AdditionalTermsEDfctU)DSG0.pluginJacgx.getRef(dsgVD).fPtr)(t, ds.getx().size(), ds.getx().getArray(), u.size(), u.getArray(), tmpJacgx.getArray(), ds.getz().size(), ds.getz().getArray());
    jacRhs += tmpJacgx; // JacRhs += \nabla_x g(x, u)
  }
  else
  {
    RuntimeException::selfThrow("ControlLinearAdditionalTermsED :: input u but no B nor pluginU");
  }
}
