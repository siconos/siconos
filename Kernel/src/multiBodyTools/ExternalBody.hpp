/* Siconos-sample version 3.1.0, Copyright INRIA 2005-2009.
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
 *
 */

/*! \file ExternalBody.hpp
  \brief External visitable lagrangian class.
*/
#ifndef ExternalBody_hpp
#define ExternalBody_hpp

#include <LagrangianDS.hpp>

DEFINE_SPTR(SpaceFilter);

class ExternalBody :
  public LagrangianDS,
  public boost::enable_shared_from_this<ExternalBody>
{
public:

  virtual void selfHash(SP::SpaceFilter) = 0;

  virtual void selfFindInteractions(SP::SpaceFilter) = 0;

  ACCEPT_VISITORS();
};

TYPEDEF_SPTR(ExternalBody);

#endif
