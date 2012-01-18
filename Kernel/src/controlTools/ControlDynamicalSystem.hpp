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

/*! \file ControlDynamicalSystem.hpp
  \brief Abstract class - General interface for all Control Dynamical Systems.
*/

#ifndef CONTROLDYNAMICALSYSTEM_H
#define CONTROLDYNAMICALSYSTEM_H

#include "SiconosPointers.hpp"
#include "TimeDiscretisation.hpp"
#include "ModelingTools.hpp"
#include "SimulationTools.hpp"

class DynamicalSystem;

class ControlDynamicalSystem : public boost::enable_shared_from_this<ControlDynamicalSystem>
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(ControlDynamicalSystem);
  /** default constructor */
  ControlDynamicalSystem() {};

protected:
  /** Constructor with the minimal set of data
   * \param a double, the starting time \f$t_0\f$
   * \param a double, the end time T
   * \param a double, the simulation time step */
  ControlDynamicalSystem(double, double, double);
  /** destructor */
  ~ControlDynamicalSystem() {};

  /** Starting time */
  double _t0;
  /** End time */
  double _T;
  /** Simulation step */
  double _h;
  /** theta for Moreau */
  double _theta;
  /** DynamicalSystem */
  SP::DynamicalSystem _processDS;
  /** Model */
  SP::Model _model;
  /** TimeDiscretisation */
  SP::TimeDiscretisation _processTD;
  /** TimeStepping */
  SP::TimeStepping _processSimulation;
  /** Moreau */
  SP::Moreau _processIntegrator;

public:
  /** Modify the value of theta (for Moreau)
   * \param an unsigned int, the new value of theta*/
  void setTheta(unsigned int);

  /** Initialize the ControlDynamicalSystem, instantiate all objects */
  void initialize();

  /** Return the _processDS */
  SP::DynamicalSystem processDS() const
  {
    return _processDS;
  };

  /** Return the Simulation */
  SP::TimeStepping simulation() const
  {
    return _processSimulation;
  };

  /** Return the Model */
  SP::Model model() const
  {
    return _model;
  };

};

DEFINE_SPTR(ControlDynamicalSystem);
#endif // CONTROLDYNAMICALSYSTEM_H
