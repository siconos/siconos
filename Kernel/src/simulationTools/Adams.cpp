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
#include "Adams.h"
using namespace std;

// --- xml constructor ---
Adams::Adams(OneStepIntegratorXML* osixml): OneStepIntegrator(osixml), r(-1)
{
  integratorType = ADAMS_INTEGRATOR;
  if (osixml != NULL)
  {
    if ((static_cast<AdamsXML*>(integratorXml))->hasR() == true)
    {
      r = (static_cast<AdamsXML*>(integratorXml))->getR();
    }
  }
  else RuntimeException::selfThrow("Adams::Adams() - xml constructor - IntegratorXML object not exists");
}

// --- Minimum data constructor ---
Adams::Adams(TimeDiscretisation* td, DynamicalSystem* ds): OneStepIntegrator(td, ds), r(-1)
{
  integratorType = ADAMS_INTEGRATOR;
}

// --- Destructor ---
Adams::~Adams()
{}

// --- Casting for Python ---
Adams* Adams::convert(OneStepIntegrator* osi)
{
  cout << "Adams::convert (OneStepIntegrator* osi)" << endl;
  Adams* adams = dynamic_cast<Adams*>(osi);
  return adams;
}

// --- Default constructor ---
Adams::Adams(): OneStepIntegrator(), r(-1)
{
  integratorType = ADAMS_INTEGRATOR;
}
