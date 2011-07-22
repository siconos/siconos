/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "TimeDiscretisation.hpp"
#include "RuntimeException.hpp"
#include "Tools.hpp"
#include <math.h>
using namespace std;

// IO Constructors -> XML
TimeDiscretisation::TimeDiscretisation(SP::TimeDiscretisationXML tdXML, double t0, double T):
  h(0.0), k(0), _timeDiscretisationXML(tdXML), tdCase(0), pos(0)
{
  if (! _timeDiscretisationXML)
    RuntimeException::selfThrow("TimeDiscretisation: xml constructor - TimeDiscretisationXML = NULL");

  // XML inputs: 3 possibilities
  //  - vector input for tk
  //  - input for h

  // --- Check what are the given data ---
  bool hasNSteps = _timeDiscretisationXML->hasN();
  bool hasH = _timeDiscretisationXML->hasH();
  bool hasTk = _timeDiscretisationXML->hasTk();

  // Eliminate cases with too many inputs
  if ((hasTk && hasH) || (hasTk && hasNSteps) || (hasH && hasNSteps))
    RuntimeException::selfThrow("TimeDiscretisation: xml constructor - Too many input data, some of them are useless.");

  // --- Read the data ---
  if (hasH) // T is useless
  {
    h = _timeDiscretisationXML->geth();
    tk.reserve(2);
    tk.push_back(t0);
    tk.push_back(t0 + h);
    tdCase = 2;
  }
  else if (hasNSteps) // t0 and T are required
  {
    unsigned int nSteps = _timeDiscretisationXML->getN();
    assert(T > t0 && "TimeDiscretisation xml constructor error: final time is less or equal to initial time.");
    h = (T - t0) / nSteps;
    tk.reserve(2);
    tk.push_back(t0);
    tk.push_back(t0 + h);
    tdCase = 2;
  }
  else if (hasTk) // neither t0 nor T is required.
  {
    // Read tk
    _timeDiscretisationXML->getTk(tk);
    h = tk[1] - tk[0];
    pos = k;
    tdCase = 1;
  }
  else
    RuntimeException::selfThrow("TimeDiscretisation: xml constructor - No input data.");

}

// --- Straightforward constructors ---

// INPUT = tk
// => tdCase = 1.
// hk is deduced from tk.
// In this case, the complete vector tk is saved in the class.
TimeDiscretisation::TimeDiscretisation(const TkVector& newTk):
  h(0.0), k(0), tdCase(1), pos(k)
{
  tk = newTk;
  h = tk[1] - tk[0];
}

// INPUTS: nSteps, t0, T  -
// => tdCase = 2. Only two values are saved for tk.
TimeDiscretisation::TimeDiscretisation(unsigned int nSteps, double t0, double T):
  h(0.0), k(0), tdCase(2), pos(0)
{
  h = (T - t0) / nSteps;
  tk.reserve(2);
  tk.push_back(t0);
  tk.push_back(t0 + h);
}

// INPUTS: t0 and h
// => tdCase = 2. Only two values are saved for tk.
TimeDiscretisation::TimeDiscretisation(double t0, double newH):
  h(newH), k(0), tdCase(2), pos(0)
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
}

void TimeDiscretisation::setTk(const TkVector& newValue)
{
  tk.clear();
  tk = newValue;
  tdCase = 1;
  pos = 0;
  k = 0;
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


  if (_timeDiscretisationXML)
  {
    _timeDiscretisationXML->setH(h);
    //_timeDiscretisationXML->setTkNode(*tk);
  }
  else RuntimeException::selfThrow("TimeDiscretisation::saveTimeDiscretisationToXML - TimeDiscretisationXML object not exists");
}

