/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

#include "NLinearBC.h"
using namespace std;

NLinearBC::NLinearBC(): BoundaryCondition()
{
  this->boundaryType = NLINEARBC;
}

NLinearBC::NLinearBC(BoundaryConditionXML* bcxml): BoundaryCondition(bcxml)
{
  this->boundaryType = NLINEARBC;
}

NLinearBC::~NLinearBC()
{}

void NLinearBC::fillBCWithBCXML()
{
  if (this->bcXML != NULL)
  {

  }
  else RuntimeException::selfThrow("NLinearBC::fillBCWithBCXML - The BoundaryConditionXML object doesn't exists");
}

void NLinearBC::saveBCToXML()
{
  if (this->bcXML != NULL)
  {

  }
  else RuntimeException::selfThrow("NLinearBC::saveBCToXML - The BoundaryConditionXML object doesn't exists");
}

void NLinearBC::createBoundaryCondition(BoundaryConditionXML * bcXML)
{
  if (bcXML != NULL)
  {
    this->bcXML = bcXML;
    this->boundaryType = NLINEARBC;
    this->fillBCWithBCXML();
  }
  else
  {}
}

NLinearBC* NLinearBC::convert(BoundaryCondition* bc)
{
  NLinearBC* nlbc = dynamic_cast<NLinearBC*>(bc);
  return nlbc;
}

