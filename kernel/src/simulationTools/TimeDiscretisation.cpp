/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "TimeDiscretisation.hpp"
#include "RuntimeException.hpp"
#include "Tools.hpp"

#include <cmath>
#include <limits>


TimeDiscretisation::TimeDiscretisation(): _h(0.), _t0(std::numeric_limits<double>::quiet_NaN())
{
  mpf_init(_hgmp);
  mpf_init(_tkp1);
  mpf_init(_tk);
  mpf_init(_t0gmp);
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
TimeDiscretisation::TimeDiscretisation(double t0, const std::string& str): _h(0.0), _t0(t0)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1);
  mpf_init(_tk);
  mpf_set_str(_hgmp, str.c_str(), 10);
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

  if(td.hGmp())
  {
    mpf_init_set(_hgmp, *td.currentTimeStep());
    _h = 0.;
  }
  else if(td.hConst())
  {
    _h = td._h;
  }
  else
  {
    _h = 0.;
  }
  _t0 = td.getT0();
  _tkV = td.getTkVector();
}


// --- Destructor ---
TimeDiscretisation::~TimeDiscretisation()
{
  if(!_tkV.empty())
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
  if(_h == 0.0)
    mpf_set_d(_t0gmp, val);
  if(!_tkV.empty())
    RuntimeException::selfThrow("TimeDiscretisation::setT0 must be called only when the TimeDiscretisation is with a constant h");
}

double TimeDiscretisation::currentTimeStep(const unsigned int k)
{
  if(_tkV.empty())
  {
    if(_h > 0.)
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
    if(_h > 0.)
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
