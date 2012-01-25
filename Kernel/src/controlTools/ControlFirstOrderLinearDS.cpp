/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

#include "ControlFirstOrderLinearDS.hpp"

using namespace std;

ControlFirstOrderLinearDS::ControlFirstOrderLinearDS(double t0, double T, double h,
    SP::SiconosVector x0, SP::SiconosMatrix A):
  ControlDynamicalSystem(t0, T, h), _x0(x0), _A(A)
{
  _processDS.reset(new FirstOrderLinearDS(_x0, _A));
}

void ControlFirstOrderLinearDS::initialize()
{
  ControlDynamicalSystem::initialize(_x0);
}
