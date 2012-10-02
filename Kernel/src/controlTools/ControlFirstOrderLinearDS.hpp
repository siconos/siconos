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

/*! \file ControlFirstOrderLinearDS.hpp
  \brief Define a FirstOrderLinearDS and set up common simulation objects in control.
*/

#ifndef CONTROLFIRSTORDERLINEARDS_H
#define CONTROLFIRSTORDERLINEARDS_H

#include "SiconosKernel.hpp"
#include "ControlDynamicalSystem.hpp"

class ControlFirstOrderLinearDS : public ControlDynamicalSystem
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(ControlFirstOrderLinearDS);

protected:
  /** Initial state \f$x_0\f$ */
  SP::SiconosVector _x0;
  /** A matrix */
  SP::SiconosMatrix _A;

  /** default constructor */
  ControlFirstOrderLinearDS() {};

public:
  /** Constructor with the minimal set of data
   * \param t0 the starting time (\f$t_0\f$)
   * \param T the end time T
   * \param h the simulation time step
   * \param x0 a SP::SiconosVector
   * \param A a SP::SiconosMatrix
   */
  ControlFirstOrderLinearDS(double t0, double T, double h, SP::SiconosVector x0, SP::SiconosMatrix A);

  /** destructor */
  virtual ~ControlFirstOrderLinearDS() {};

  /** Initialization */
  void initialize();

  /** Return the _processDS */
  SP::FirstOrderLinearDS processDS() const
  {
    return std11::static_pointer_cast<FirstOrderLinearDS>(_processDS);
  };
};

DEFINE_SPTR(ControlFirstOrderLinearDS)
#endif // CONTROLFIRSTORDERLINEARDS_H
