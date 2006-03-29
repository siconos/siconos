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
#include "OneStepIntegrator.h"
using namespace std;

// --- Xml constructor ---
OneStepIntegrator::OneStepIntegrator(OneStepIntegratorXML* osixml): integratorType("undefined"), ds(NULL), sizeMem(1), timeDiscretisation(NULL), integratorXml(osixml)
{
  if (integratorXml != NULL)
  {
    if (integratorXml->hasR()) sizeMem = integratorXml->getR();
  }
  else RuntimeException::selfThrow("OneStepIntegrator: xml constructor, xml object = NULL");
}

// --- Constructor from a minimum set of data ---
OneStepIntegrator::OneStepIntegrator(TimeDiscretisation* td, DynamicalSystem* newDs): integratorType("undefined"), ds(newDs), sizeMem(1), timeDiscretisation(td), integratorXml(NULL)
{}

// --- Destructor ---
OneStepIntegrator::~OneStepIntegrator()
{
  timeDiscretisation = NULL;
  ds = NULL;
  integratorXml = NULL;
}


void OneStepIntegrator::initialize()
{
  double t0 = timeDiscretisation->getT0();
  ds->initialize(t0, sizeMem);
}

void OneStepIntegrator::nextStep()
{
  ds->setIsDSUp(false); // to reset isDSUp bool, see LagrangianDS.
  ds->swapInMemory();
  ds->getRPtr()->zero();
}

void OneStepIntegrator::computeFreeState()
{
  RuntimeException::selfThrow("OneStepIntegrator:computeFreeState, not yet implemented for this type of integrator" + getType());
}

void OneStepIntegrator::integrate()
{
  RuntimeException::selfThrow("OneStepIntegrator:integrate, not yet implemented for this type of integrator" + getType());
}

void OneStepIntegrator::integrate(const double&, const double&, const double&, const bool&)
{

  cout << " OSI Integrate, not yet implemented, test version" << endl;
}

void OneStepIntegrator::updateState()
{
  RuntimeException::selfThrow("OneStepIntegrator:updateState, not yet implemented for this type of integrator" + getType());
}

void OneStepIntegrator::display() const
{
  cout << "==== OneStepIntegrator display =====" << endl;
  cout << "| integratorType : " << integratorType << endl;
  cout << "| DS id is: " << endl;
  if (ds != NULL) cout << ds->getId() << endl;
  else cout << "-> NULL" << endl;
  cout << "| sizeMem: " << sizeMem << endl;
  if (timeDiscretisation != NULL) timeDiscretisation->display();
  else cout << "-> NULL" << endl;
  cout << "====================================" << endl;
}

void OneStepIntegrator::saveIntegratorToXML()
{
  IN("OneStepIntegrator::saveIntegratorToXML\n");
  if (integratorXml != 0)
  {
    vector<int> dsConcerned;
    dsConcerned.push_back(ds->getNumber());
    integratorXml->setDSConcerned(&dsConcerned);

    // r is saved only if the integrator is not a Moreau integrator !
    if (integratorType != MOREAU_INTEGRATOR) integratorXml->setR(sizeMem);
  }
  else RuntimeException::selfThrow("OneStepIntegrator::saveIntegratorToXML - OneStepIntegratorXML object = NULL");
  OUT("OneStepIntegrator::saveIntegratorToXML\n");

}

//-- Default constructor --
OneStepIntegrator::OneStepIntegrator(): integratorType("none"), ds(NULL), sizeMem(1), timeDiscretisation(NULL), integratorXml(NULL)
{}
