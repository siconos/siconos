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
  _h(0.0), _k(0), _timeDiscretisationXML(tdXML), _tdCase(0), _pos(0)
{
  if (!_timeDiscretisationXML)
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
    _h = _timeDiscretisationXML->geth();
    _tk.reserve(2);
    _tk.push_back(t0);
    _tk.push_back(t0 + _h);
    _tdCase = 2;
  }
  else if (hasNSteps) // t0 and T are required
  {
    unsigned int nSteps = _timeDiscretisationXML->getN();
    assert(T > t0 && "TimeDiscretisation xml constructor error: final time is less or equal to initial time.");
    _h = (T - t0) / nSteps;
    _tk.reserve(2);
    _tk.push_back(t0);
    _tk.push_back(t0 + _h);
    _tdCase = 2;
  }
  else if (hasTk) // neither t0 nor T is required.
  {
    // Read tk
    _timeDiscretisationXML->getTk(_tk);
    _h = _tk[1] - _tk[0];
    _pos = _k;
    _tdCase = 1;
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
  _h(0.0), _k(0), _tdCase(1), _pos(_k)
{
  _tk = newTk;
  _h = _tk[1] - _tk[0];
}

// INPUTS: nSteps, t0, T
// => tdCase = 2. Only two values are saved for tk.
TimeDiscretisation::TimeDiscretisation(unsigned int nSteps, double t0, double T):
  _h(0.0), _k(0), _tdCase(2), _pos(0)
{
  _h = (T - t0) / nSteps;
  _tk.reserve(2);
  _tk.push_back(t0);
  _tk.push_back(t0 + _h);
}

// INPUTS: t0 and h
// => tdCase = 2. Only two values are saved for tk.
TimeDiscretisation::TimeDiscretisation(double t0, double newH):
  _h(newH), _k(0), _tdCase(2), _pos(0)
{
  _tk.reserve(2);
  _tk.push_back(t0);
  _tk.push_back(t0 + _h);
}


// Copy constructor
TimeDiscretisation::TimeDiscretisation(const TimeDiscretisation& td)
{
  _h = td.currentTimeStep();
  _k = td.getK();
  _tdCase = td.getTDCase();
  // the magic formula is in the header file
  _pos = _tdCase == 1 ? _k : 0;
  _tk = td.getTk();
  if (td.timeDiscretisationXML())
    _timeDiscretisationXML.reset(new TimeDiscretisationXML(*(td.timeDiscretisationXML())));
}


// --- Destructor ---
TimeDiscretisation::~TimeDiscretisation()
{
  _tk.clear();
}

void TimeDiscretisation::setCurrentTimeStep(double newH)
{
  _h = newH;
  if (_tdCase == 1)
  {
    cout << "TimeDiscretisation::setCurrentTimeStep(newH) Warning: you change the time step whereas the Time Discretisation was built with a complete tk vector. This will result in a change in all future time steps: from now on, h = newH until next call to this function. " << endl;
    _tdCase = 2;
    _tk[0] = _tk[_k];
    _pos = 0;
    _tk.resize(2);
  }
  _tk[_pos + 1] = _tk[_pos] + _h;
}

void TimeDiscretisation::setTk(const TkVector& newTk)
{
  _tk.clear();
  _tk = newTk;
  _tdCase = 1;
  _pos = 0;
  _k = 0;
  _h = _tk[1] - _tk[0];
}

void TimeDiscretisation::increment()
{
  _k++;
  if (_tdCase == 1) // h is deduced from tk
  {
    _pos = _k;
    _h = _tk[_pos + 1] - _tk[_pos];
  }
  else // tk+1 is build with tk and h
  {
    _tk[_pos] = _tk[_pos + 1];
    _tk[_pos + 1] = _tk[_pos] + _h;
  }
}

// --- Other functions ---
void TimeDiscretisation::display() const
{
  cout << "====> Time Disretisation :" << endl;
  cout << " current time step starts from " << _tk[_pos] << " and ends at " << _tk[_pos + 1] << endl;
  cout << "====" << endl;
}

// --- XML functions ---

void TimeDiscretisation::saveTimeDiscretisationToXML()
{
  RuntimeException::selfThrow("TimeDiscretisation::saveTimeDiscretisationToXML -Not yet properly implemented");

  if (_timeDiscretisationXML)
  {
    _timeDiscretisationXML->setH(_h);
    //_timeDiscretisationXML->setTkNode(*tk);
  }
  else RuntimeException::selfThrow("TimeDiscretisation::saveTimeDiscretisationToXML - TimeDiscretisationXML object not exists");
}

