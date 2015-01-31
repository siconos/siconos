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

/*!\file ControlLinearAdditionalTermsTS.hpp
 * \brief Functions to add control terms during the integration steps with a TimeStepping scheme
 */

#include "ExtraAdditionalTerms.hpp"

struct ControlLinearAdditionalTermsTS : ExtraAdditionalTerms
{

  /** initialize elements in the graph for the computations
   * \param DSG0 the graph of DynamicalSystems
   * \param model the current Model
   */
  virtual void init(DynamicalSystemsGraph& DSG0, const Model& model);

  /** add smooth term to xfree (like the control input, the error correction for an observer)
   * \param DSG0 the graph of DynamicalSystems
   * \param dsgVD a DynamicalSystem in the DS graph
   * \param h the current timestep
   * \param xfree the free state to modify
   */
  virtual void addSmoothTerms(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double h, SiconosVector& xfree);

  /** add contribution to JacRhs for instance if \f$\dot{x} = f(x) + g(x)u\f$
   * \param DSG0 the graph of DynamicalSystems
   * \param dsgVD a DynamicalSystem in the DS graph
   * \param t the current timestep
   * \param jacRhs the jacobian to modify
   */
  virtual void addJacobianRhsContribution(DynamicalSystemsGraph& DSG0, const DynamicalSystemsGraph::VDescriptor& dsgVD, const double t, SiconosMatrix& jacRhs);

};
