/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP


#include <vector>
#include "SiconosPointers.hpp"
#include "SimpleVector.hpp"
#include "PluggedObject.hpp"
typedef  void (*FPtrPrescribedVelocity)(double, unsigned int, double*);

using namespace std;

/** \class BoundaryCondition
 *  \brief This class models simple boundary conditions for prescribing the velocities
 *   in a Dynamical System. A simple boundary condition is considered to fix a component \f$ j \f$ of *   the velocity vector, i.e., \f$ v_j(t) = bc(t)\f$ where \f$ bc(t)\f$ is a given function of time.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.2.0.
 *  \date (Re-Creation) September 2010
 *
 *
 */
class BoundaryCondition
{
public:

  /** \fn BoundaryCondition();
   *  \brief Basic constructor
   */
  //   BoundaryCondition();

  BoundaryCondition(std::vector<unsigned int> * newVelocityIndices,  SP::SimpleVector newVelocityValues);

  virtual ~BoundaryCondition();

  inline vector<unsigned int>  * velocityIndices()
  {
    return _velocityIndices;
  };

  inline SP::SimpleVector prescribedVelocity()
  {
    return _prescribedVelocity;
  };

  inline SP::SimpleVector reactionImpulse()
  {
    return _reactionImpulse;
  };

  /** allow to set a specified function to compute prescribedVelocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputePrescribedVelocityFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginPrescribedVelocity->setComputeFunction(pluginPath, functionName);
    if (!_prescribedVelocity) _prescribedVelocity.reset(new SimpleVector(_velocityIndices->size()));
  }

  /** default function to compute the external strengths
   *  \param double time : the current time
   */
  virtual void computePrescribedVelocity(double time);

protected:
  /* Indices of the prescribed component of the velocity vector */
  std::vector<unsigned int>   * _velocityIndices;
  /* Values of the prescribed component of the velocity vector */
  SP::SimpleVector _prescribedVelocity;
  /* Values of the multiplier (impulse) associted with */
  SP::SimpleVector _reactionImpulse;
  /*plugin defining the function V(t)*/
  SP::PluggedObject _pluginPrescribedVelocity;
  //   /*Link to the precribed DynamicalSystem*/
  //   SP::DynamicalSystem _DS;
};
TYPEDEF_SPTR(BoundaryCondition);
#endif // BOUNDARYCONDITION_HPP
