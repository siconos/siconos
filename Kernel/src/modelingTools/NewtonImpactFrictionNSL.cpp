/* Siconos version 1.0, Copyright INRIA 2005.
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
#include "NewtonImpactFrictionNSL.h"
using namespace std;

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL():
  NonSmoothLaw(), en(0.0), et(0.0), mu(0.0)
{
  nsLawType = NEWTONIMPACTFRICTIONNSLAW;
}

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(NonSmoothLawXML* nslawxml):
  NonSmoothLaw(nslawxml), en(0.0), et(0.0), mu(0.0)
{
  nsLawType = NEWTONIMPACTFRICTIONNSLAW;
  if (nslawxml != NULL)
  {
    en = (static_cast<NewtonImpactFrictionNSLXML*>(nslawxml))->getEn();
    et = (static_cast<NewtonImpactFrictionNSLXML*>(nslawxml))->getEt();
    mu = (static_cast<NewtonImpactFrictionNSLXML*>(nslawxml))->getMu();
  }
  else RuntimeException::selfThrow("NewtonImpactFrictionNSL:: xml constructor, xml file=NULL");
}

NewtonImpactFrictionNSL::NewtonImpactFrictionNSL(const double& newEn, const double& newEt, const double& newMu):
  NonSmoothLaw(), en(newEn), et(newEt), mu(newMu)
{
  nsLawType = NEWTONIMPACTFRICTIONNSLAW;
}

NewtonImpactFrictionNSL::~NewtonImpactFrictionNSL()
{}

bool NewtonImpactFrictionNSL::isVerified() const
{
  bool res = false;
  // to do
  RuntimeException::selfThrow("NewtonImpactFrictionNSL:: isVerified, not yet implemented!");
  return res;
}

void NewtonImpactFrictionNSL::display() const
{
  cout << "=== Newton impact-friction non-smooth law data display ===" << endl;
  cout << " Normal Newton coefficient of restitution: " << en << endl;
  cout << " Tangential Newton coefficient of restitution: " << et << endl;
  cout << "Friction coefficient: " << mu << endl;
  cout << "==========================================================" << endl;
}

void NewtonImpactFrictionNSL::saveNonSmoothLawToXML()
{
  IN("NewtonImpactFrictionNSL::saveNonSmoothLawToXML\n");
  static_cast<NewtonImpactFrictionNSLXML*>(nslawxml)->setEn(en);
  static_cast<NewtonImpactFrictionNSLXML*>(nslawxml)->setEt(et);
  static_cast<NewtonImpactFrictionNSLXML*>(nslawxml)->setMu(mu);
  OUT("NewtonImpactFrictionNSL::saveNonSmoothLawToXML\n");
}

NewtonImpactFrictionNSL* NewtonImpactFrictionNSL::convert(NonSmoothLaw* nsl)
{
  cout << "NewtonImpactFrictionNSL::convert (NonSmoothLaw* nsl)" << endl;
  NewtonImpactFrictionNSL* nilnsl = dynamic_cast<NewtonImpactFrictionNSL*>(nsl);
  return nilnsl;
}


