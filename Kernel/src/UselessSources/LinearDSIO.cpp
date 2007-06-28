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
#include "LinearDSIO.h"
using namespace std;

LinearDSIO::LinearDSIO(): DSInputOutput()
{
  this->dsioType = LINEARDSIO;
}

LinearDSIO::LinearDSIO(DSInputOutputXML* dsioxml): DSInputOutput(dsioxml)
{
  this->dsioType = LINEARDSIO;
}

LinearDSIO::~LinearDSIO()
{
  if (A != NULL) delete A;
  if (B != NULL) delete B;
}


void LinearDSIO::fillDSInputOutputWithDSInputOutputXML()
{
  if (this->dsioxml != NULL)
  {
    this->number = this->dsioxml->getNumber();
    A = new SimpleMatrix(static_cast<LinearDSIOXML*>(this->dsioxml)->getA());
    B = new SimpleMatrix(static_cast<LinearDSIOXML*>(this->dsioxml)->getB());
  }
  else RuntimeException::selfThrow("DSInputOutput::fillDSInputOutputWithDSInputOutputXML - object DSInputOutputXML does not exist");
}

void LinearDSIO::saveDSInputOutputToXML()
{
  if (this->dsioxml != NULL)
  {
    /*
     * these attributes are only required for LagrangianNonLinear DSInputOutput !
     */
    static_cast<LinearDSIOXML*>(this->dsioxml)->setA(A);
    static_cast<LinearDSIOXML*>(this->dsioxml)->setB(B);
  }
  else RuntimeException::selfThrow("DSInputOutput::saveDSInputOutputToXML - object DSInputOutputXML does not exist");
}

void LinearDSIO::createDSInputOutput(DSInputOutputXML * dsioXML, int number,
                                     SiconosMatrix *newA, SiconosMatrix *newB)
{
  if (dsioXML != NULL)
  {
    //    this->init();
    this->dsioxml = dsioXML;
    this->dsioType = LINEARDSIO;
    this->fillDSInputOutputWithDSInputOutputXML();
  }
  else
  {
    //this->dsioxml = dsioXML;
    this->dsioType = LINEARDSIO;
    this->number = number;
    A = new SimpleMatrix(*A);
    B = new SimpleMatrix(*B);
  }
}

