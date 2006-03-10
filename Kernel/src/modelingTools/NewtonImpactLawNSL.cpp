/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#include "NewtonImpactLawNSL.h"
using namespace std;

NewtonImpactLawNSL::NewtonImpactLawNSL(): NonSmoothLaw(), e(0.0)
{
  nsLawType = NEWTONIMPACTNSLAW;
}

NewtonImpactLawNSL::NewtonImpactLawNSL(NonSmoothLawXML* nslawxml):
  NonSmoothLaw(nslawxml), e(0.0)
{
  nsLawType = NEWTONIMPACTNSLAW;
  if (nslawxml != NULL)
    e = (static_cast<NewtonImpactLawNSLXML*>(nslawxml))->getE();
  else RuntimeException::selfThrow("NewtonImpactLawNSL:: xml constructor, xml file=NULL");
}

NewtonImpactLawNSL::NewtonImpactLawNSL(const double& newE):
  NonSmoothLaw(), e(newE)
{
  nsLawType = NEWTONIMPACTNSLAW;
}

NewtonImpactLawNSL::~NewtonImpactLawNSL()
{}

bool NewtonImpactLawNSL::isVerified() const
{
  bool res = false;
  // to do
  RuntimeException::selfThrow("NewtonImpactFrictionNSL:: isVerified, not yet implemented!");
  return res;
}

void NewtonImpactLawNSL::display() const
{
  cout << "===============================================================================" << endl;
  cout << "=== Newton impact (frictionless) non-smooth law coefficient of restitution: " << e << endl;
  cout << "===============================================================================" << endl;
}

void NewtonImpactLawNSL::saveNonSmoothLawToXML()
{
  IN("NewtonImpactLawNSL::saveNonSmoothLawToXML\n");
  static_cast<NewtonImpactLawNSLXML*>(this->nslawxml)->setE(e);
  OUT("NewtonImpactLawNSL::saveNonSmoothLawToXML\n");
}

NewtonImpactLawNSL* NewtonImpactLawNSL::convert(NonSmoothLaw* nsl)
{
  NewtonImpactLawNSL* nilnsl = dynamic_cast<NewtonImpactLawNSL*>(nsl);
  return nilnsl;
}


