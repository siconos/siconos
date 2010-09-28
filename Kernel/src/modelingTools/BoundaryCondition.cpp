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

#include "BoundaryCondition.hpp"

using namespace std;

// BoundaryCondition::BoundaryCondition()
// _velocityIndices(NULL)
// {

// }

BoundaryCondition::BoundaryCondition(std::vector<unsigned int> * newVelocityIndices, SP::SimpleVector newVelocityValues): _velocityIndices(newVelocityIndices),  _prescribedVelocity(newVelocityValues)
{


}

// BoundaryCondition::initialize(SP::DynamicalSystem newDS)
// {
//   _DS = newDS;
//   int numberofBoundaryConditions = _velocityIndices.size();

//   if (_DS->getType() == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
//   {
//     _Wcolumns.reset(new SimpleMatrix(_DS->getNdof ,numberofBoundaryConditions ))  ;
//   }
//   else
//   {
//     _Wcolumns.reset(new SimpleMatrix(_DS->getN() , numberofBoundaryConditions ))  ;
//   }

// //   vector<unsigned int>::iterator itindex;
// //   for (itindex =_velocityIndices.begin() ;itindex =_velocityIndices.end(); itindex++ )
// //   {
// //   }

// }


BoundaryCondition::~BoundaryCondition()
{
}

void BoundaryCondition::computePrescribedVelocity(double time)
{
  if (_pluginPrescribedVelocity->fPtr)
    ((FPtrPrescribedVelocity)_pluginPrescribedVelocity->fPtr)(time, _velocityIndices->size(), &(*_prescribedVelocity)(0));
}
