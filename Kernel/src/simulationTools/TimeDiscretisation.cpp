/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "TimeDiscretisation.h"
#include "RuntimeException.h"
#include "Model.h"
#include "Tools.h"
#include <math.h>
using namespace std;

// IO Constructors -> XML
TimeDiscretisation::TimeDiscretisation(SP::TimeDiscretisationXML tdXML, SP::Model m):
  h(0.0), k(0), timeDiscretisationXML(tdXML), model(m), isUpToDate(false), tdCase(0), pos(0)
{
  if (! timeDiscretisationXML)
    RuntimeException::selfThrow("TimeDiscretisation: xml constructor - TimeDiscretisationXML = NULL");

  // XML inputs: 3 possibilities
  //  - vector input for tk
  //  - input for h

  // --- Check what are the given data ---
  bool hasNSteps = timeDiscretisationXML->hasN();
  bool hasH = timeDiscretisationXML->hasH();
  bool hasTk = timeDiscretisationXML->hasTk();

  // Eliminate cases with too many inputs
  if ((hasTk && hasH) || (hasTk && hasNSteps) || (hasH && hasNSteps))
    RuntimeException::selfThrow("TimeDiscretisation: xml constructor - Too many input data, some of them are useless.");

  // --- Read the data ---
  if (hasH)
  {
    h = timeDiscretisationXML->getH();
    if (!model)
      RuntimeException::selfThrow("TimeDiscretisation xml constructor - Not linked to a valid Model (NULL pointer).");

    double t0 = model->getT0();
    double T = model->getFinalT();
    tk.reserve(2);
    tk.push_back(t0);
    tk.push_back(t0 + h);
    isUpToDate = true;
    tdCase = 2;
  }
  else if (hasNSteps)
  {
    unsigned int nSteps = timeDiscretisationXML->getN();
    if (!model)
      RuntimeException::selfThrow("TimeDiscretisation xml constructor - Not linked to a valid Model (NULL pointer).");

    double t0 = model->getT0();
    double T = model->getFinalT();
    h = (T - t0) / nSteps;
    tk.reserve(2);
    tk.push_back(t0);
    tk.push_back(t0 + h);
    isUpToDate = true;
    tdCase = 2;
  }
  else if (hasTk)
  {
    // Read tk
    timeDiscretisationXML->getTk(tk);
    h = tk[1] - tk[0];
    pos = k;
    tdCase = 1;
  }
  else
    RuntimeException::selfThrow("TimeDiscretisation: xml constructor - No input data.");

}

// --- Straightforward constructors ---

// INPUT = tk - Model is optional.
// => tdCase = 1.
// hk is deduced from tk.
// In this case, the complete vector tk is saved in the class.
TimeDiscretisation::TimeDiscretisation(const TkVector& newTk, SP::Model m):
  h(0.0), k(0), model(m), isUpToDate(true), tdCase(1), pos(k)
{
  tk = newTk;
  h = tk[1] - tk[0];
}

// INPUTS: nSteps - Model is required (for t0 and T) -
// => tdCase = 2. Only two values are saved for tk.
TimeDiscretisation::TimeDiscretisation(unsigned int nSteps, SP::Model m):
  h(0.0), k(0), model(m), isUpToDate(true), tdCase(2), pos(0)
{
  if (!model)
    RuntimeException::selfThrow("TimeDiscretisation::initialize - Not linked to a valid Model (NULL pointer).");

  double t0 = model->getT0();
  double T = model->getFinalT();
  h = (T - t0) / nSteps;
  tk.reserve(2);
  tk.push_back(t0);
  tk.push_back(t0 + h);
}

// INPUTS: h - Model is required (for t0 and T) -
// => tdCase = 2. Only two values are saved for tk.
TimeDiscretisation::TimeDiscretisation(double newH, SP::Model m):
  h(newH), k(0), model(m), isUpToDate(true), tdCase(2), pos(0)
{
  if (!model)
    RuntimeException::selfThrow("TimeDiscretisation constructor - Not linked to a valid Model (NULL pointer).");

  double t0 = model->getT0();
  tk.reserve(2);
  tk.push_back(t0);
  tk.push_back(t0 + h);
}

// INPUTS: t0 and h
// => tdCase = 2. Only two values are saved for tk.
TimeDiscretisation::TimeDiscretisation(double t0, double newH):
  h(newH), k(0), isUpToDate(true), tdCase(2), pos(0)
{
  tk.reserve(2);
  tk.push_back(t0);
  tk.push_back(t0 + h);
}

// --- Destructor ---
TimeDiscretisation::~TimeDiscretisation()
{
  tk.clear();
}

void TimeDiscretisation::setCurrentTimeStep(double newH)
{
  h = newH;
  if (tdCase == 1)
  {
    cout << "TimeDiscretisation::setCurrentTimeStep(newH) Warning: you change the time step whereas the Time Discretisation was built with a complete tk vector. This will result in a change in all future time steps: from now on, h = newH until next call to this function. " << endl;
    tdCase = 2;
    tk[0] = tk[k];
    pos = 0;
    tk.resize(2);
  }
  tk[pos + 1] = tk[pos] + h;
  isUpToDate = true;
}

void TimeDiscretisation::setTk(const TkVector& newValue)
{
  tk.clear();
  tk = newValue;
  tdCase = 1;
  pos = 0;
  k = 0;
  isUpToDate = true;
  h = tk[1] - tk[0];
}

void TimeDiscretisation::increment()
{
  k++;
  if (tdCase == 1) // h is deduced from tk
  {
    pos = k;
    h = tk[pos + 1] - tk[pos];
  }
  else // tk+1 is build with tk and h
  {
    tk[pos] = tk[pos + 1];
    tk[pos + 1] = tk[pos] + h;
  }
}

// --- Other functions ---
void TimeDiscretisation::display() const
{
  cout << "====> Time Disretisation :" << endl;
  cout << " current time step starts from " << tk[pos] << " and ends at " << tk[pos + 1] << endl;
  cout << "====" << endl;
}

// --- XML functions ---

void TimeDiscretisation::saveTimeDiscretisationToXML()
{
  RuntimeException::selfThrow("TimeDiscretisation::saveTimeDiscretisationToXML -Not yet properly implemented");


  if (timeDiscretisationXML)
  {
    timeDiscretisationXML->setH(h);
    //timeDiscretisationXML->setTkNode(*tk);
  }
  else RuntimeException::selfThrow("TimeDiscretisation::saveTimeDiscretisationToXML - TimeDiscretisationXML object not exists");
}

