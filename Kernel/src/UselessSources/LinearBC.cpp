/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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

#include "LinearBC.h"
using namespace std;

LinearBC::LinearBC(): BoundaryCondition()
{
  this->boundaryType = LINEARBC;
}

LinearBC::LinearBC(BoundaryConditionXML* bcxml): BoundaryCondition(bcxml)
{
  this->boundaryType = LINEARBC;
}

LinearBC::~LinearBC()
{
  if (omega != NULL) delete omega;
  if (omega0 != NULL) delete omega0;
  if (omegaT != NULL) delete omegaT;
}


void LinearBC::fillBCWithBCXML()
{
  if (this->bcXML != NULL)
  {
    omega = new SimpleVector(static_cast<LinearBCXML*>(this->bcXML)->getOmega());
    omegaT = new SimpleMatrix(static_cast<LinearBCXML*>(this->bcXML)->getOmegaT());
    omega0 = new SimpleMatrix(static_cast<LinearBCXML*>(this->bcXML)->getOmega0());
  }
  else RuntimeException::selfThrow("LinearBC::fillBCWithBCXML - The BoundaryConditionXML object doesn't exists");
}

void LinearBC::saveBCToXML()
{
  if (this->bcXML != NULL)
  {
    static_cast<LinearBCXML*>(this->bcXML)->setOmega(omega);
    static_cast<LinearBCXML*>(this->bcXML)->setOmegaT(omegaT);
    static_cast<LinearBCXML*>(this->bcXML)->setOmega0(omega0);
  }
  else RuntimeException::selfThrow("LinearBC::saveBCToXML - The BoundaryConditionXML object doesn't exists");
}

void LinearBC::createBoundaryCondition(BoundaryConditionXML * bcXML,
                                       SiconosVector* newOmega, SiconosMatrix* newOmega0, SiconosMatrix* newOmegaT)
{
  if (bcXML != NULL)
  {
    this->bcXML = bcXML;
    this->boundaryType = LINEARBC;
    this->fillBCWithBCXML();
  }
  else if (omega != NULL && omega0 != NULL && omegaT != NULL)
  {
    omega = new SimpleVector(*newOmega);
    omega0 = new SimpleMatrix(*newOmega0);
    omegaT = new SimpleMatrix(*newOmegaT);
  }
  else RuntimeException::selfThrow("LinearBC::createBoundaryCondition - The omega, omega0 and/or omegaT matrices is/are missing");
}


LinearBC* LinearBC::convert(BoundaryCondition* bc)
{
  LinearBC* lbc = dynamic_cast<LinearBC*>(bc);
  return lbc;
}

