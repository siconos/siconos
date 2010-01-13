/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

#include "RelayNSL.hpp"
#include "RelayNSLXML.hpp"
using namespace std;

// Default (private)
RelayNSL::RelayNSL(): NonSmoothLaw(), _lb(-1.0), _ub(1.0)
{}

RelayNSL::RelayNSL(unsigned int newSize):
  NonSmoothLaw(newSize), _lb(-1.0), _ub(1.0)
{}

RelayNSL::RelayNSL(SP::NonSmoothLawXML nslawxml):
  NonSmoothLaw(nslawxml), _lb(-1.0), _ub(1.0)
{
  if (nslawxml)
  {
    _lb = (boost::static_pointer_cast<RelayNSLXML>(nslawxml))->getC();
    _ub = (boost::static_pointer_cast<RelayNSLXML>(nslawxml))->getD();
  }
  else RuntimeException::selfThrow("RelayNSL::xml constructor, xml file=NULL");
}

RelayNSL::RelayNSL(double newLb, double newUb, unsigned int newSize):
  NonSmoothLaw(newSize), _lb(newLb), _ub(newUb)
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
  cout << "| nSLawSize : " << _size << endl;
  cout << "| lb : " << _lb << endl;
  cout << "| ub : " << _ub << endl;
  cout << "____________________________" << endl;
  cout << "------------------------------------" << endl;
}

void RelayNSL::saveNonSmoothLawToXML()
{
  boost::static_pointer_cast<RelayNSLXML>(_nslawxml)->setC(_lb);
  boost::static_pointer_cast<RelayNSLXML>(_nslawxml)->setD(_ub);
}

RelayNSL* RelayNSL::convert(NonSmoothLaw* nsl)
{
  RelayNSL* rnsl = dynamic_cast<RelayNSL*>(nsl);
  return rnsl;
}
