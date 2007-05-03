/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

#include "RelayNSL.h"
#include "RelayNSLXML.h"
using namespace std;

// Default (private)
RelayNSL::RelayNSL(): NonSmoothLaw(), c(0.0), d(0.0)
{}

RelayNSL::RelayNSL(unsigned int newSize):
  NonSmoothLaw(RELAYNSLAW, newSize), c(0.0), d(0.0)
{}

RelayNSL::RelayNSL(NonSmoothLawXML* nslawxml):
  NonSmoothLaw(RELAYNSLAW, nslawxml), c(0.0), d(0.0)
{
  if (nslawxml != NULL)
  {
    c = (static_cast<RelayNSLXML*>(nslawxml))->getC();
    d = (static_cast<RelayNSLXML*>(nslawxml))->getD();
  }
  else RuntimeException::selfThrow("RelayNSL::xml constructor, xml file=NULL");
}

RelayNSL::RelayNSL(double newC, double newD, unsigned int newSize):
  NonSmoothLaw(RELAYNSLAW, newSize), c(newC), d(newD)
{}

RelayNSL::~RelayNSL()
{}

bool RelayNSL::isVerified(void) const
{
  bool res = false;
  // to do
  return res;
}

void RelayNSL::display() const
{
  cout << "------------------------------------" << endl;
  cout << "____ data of the RelayNSL" << endl;
  cout << "| c : " << c << endl;
  cout << "| d : " << d << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void RelayNSL::saveNonSmoothLawToXML()
{
  static_cast<RelayNSLXML*>(nslawxml)->setC(c);
  static_cast<RelayNSLXML*>(nslawxml)->setD(d);
}

RelayNSL* RelayNSL::convert(NonSmoothLaw* nsl)
{
  RelayNSL* rnsl = dynamic_cast<RelayNSL*>(nsl);
  return rnsl;
}
