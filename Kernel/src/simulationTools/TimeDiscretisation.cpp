/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include <cmath>

TimeDiscretisation::TimeDiscretisation()
{
  mpf_init(_hgmp);
  mpf_init(_tkp1); 
  mpf_init(_tk); 
  mpf_init(_t0gmp);
}

// IO Constructors -> XML
TimeDiscretisation::TimeDiscretisation(SP::TimeDiscretisationXML tdXML, double t0, double T):
  _h(0.0), _timeDiscretisationXML(tdXML), _t0(t0)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1); 
  mpf_init(_tk); 
  mpf_init(_t0gmp);
  
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
  }
  else if (hasNSteps) // t0 and T are required
  {
    unsigned int nSteps = _timeDiscretisationXML->getN();
    assert(T > t0 && "TimeDiscretisation xml constructor error: final time is less or equal to initial time.");
    _h = (T - t0) / nSteps;
  }
  else if (hasTk) // neither t0 nor T is required.
  {
    // Read tk
    _timeDiscretisationXML->getTk(_tkV);
    _h = _tkV[1] - _tkV[0];
  }
  else
    RuntimeException::selfThrow("TimeDiscretisation: xml constructor - No input data.");

}

// --- Straightforward constructors ---

TimeDiscretisation::TimeDiscretisation(const TkVector& tk):
  _h(0.0)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1); 
  mpf_init(_tk); 
  mpf_init(_t0gmp);

  _tkV = tk;
  _t0 = _tkV.at(0);
}

// INPUTS: t0 and h
TimeDiscretisation::TimeDiscretisation(double t0, double h):
  _h(h), _t0(t0)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1); 
  mpf_init(_tk); 
  mpf_init(_t0gmp);

}

// INPUTS: t0 and h
TimeDiscretisation::TimeDiscretisation(double t0, std::string& str): _t0(t0)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1);
  mpf_init(_tk);
  mpf_init(_t0gmp);
  mpf_set_str(_hgmp, str.c_str(), 10);
  _h = 0.0;
  mpf_init_set_d(_t0gmp, t0);
}

TimeDiscretisation::TimeDiscretisation(unsigned int nSteps, double t0, double T):
  _t0(t0)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1); 
  mpf_init(_tk); 
  mpf_init(_t0gmp);

  _h = (T - t0) / nSteps;
}

// Copy constructor
TimeDiscretisation::TimeDiscretisation(const TimeDiscretisation& td)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1);
  mpf_init(_tk);
  mpf_init(_t0gmp);
  
  if (td.hGmp())
  {
    mpf_init_set(_hgmp, *td.currentTimeStep());
  }
  else
  {
    if (td.hConst())
      _h = td._h;
  }
  _t0 = td.getT0();
  _tkV = td.getTkVector();
  if (td.timeDiscretisationXML())
    _timeDiscretisationXML.reset(new TimeDiscretisationXML(*(td.timeDiscretisationXML())));
}


// --- Destructor ---
TimeDiscretisation::~TimeDiscretisation()
{
  if (!_tkV.empty())
    _tkV.clear();
  
  mpf_clear(_hgmp);
  mpf_clear(_tkp1);
  mpf_clear(_tk);
  mpf_clear(_t0gmp);

}

void TimeDiscretisation::setTkVector(const TkVector& newTk)
{
  _tkV.clear();
  _tkV = newTk;
}

void TimeDiscretisation::setT0(double val)
{
  _t0 = val;
  if (_h == 0.0)
    mpf_set_d(_t0gmp, val);
  if (!_tkV.empty())
    RuntimeException::selfThrow("TimeDiscretisation::setT0 must be called only when the TimeDiscretisation is with a constant h");
}

double TimeDiscretisation::currentTimeStep(const unsigned int k)
{
  if(_tkV.empty())
  {
    if (_h > 0)
      return _h;
    else
    {
      mpf_mul_ui(_tkp1, _hgmp, k+1);
      mpf_mul_ui(_tk, _hgmp, k);
      mpf_add(_tk, _tk, _t0gmp);
      mpf_add(_tkp1, _tkp1, _t0gmp);
      return mpf_get_d(_tkp1) - mpf_get_d(_tk);
    }
  }
  else
    return _tkV.at(k+1) - _tkV.at(k);
}

double TimeDiscretisation::getTk(const unsigned int indx)
{
  if(_tkV.empty())
  {
    if (_h > 0)
      return _t0 + _h*indx;
    else
    {
      mpf_mul_ui(_tk, _hgmp, indx);
      mpf_add(_tk, _tk, _t0gmp);
      return mpf_get_d(_tk);
    }
  }
  else
    return _tkV.at(indx);
}

// --- Other functions ---
void TimeDiscretisation::display() const
{
  std::cout << "====> Time Disretisation :" <<std::endl;
  std::cout << " the current timestep is " << _h << std::endl;
  std::cout << "====" <<std::endl;
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

