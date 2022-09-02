/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include <cmath>
#include <limits>

#include "SiconosException.hpp"
#include "Tools.hpp"

siconos::simulation::TimeDiscretisation::TimeDiscretisation(double t0, double h) : _h(h), _t0(t0)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1);
  mpf_init(_tk);
  mpf_init(_t0gmp);
}

siconos::simulation::TimeDiscretisation::TimeDiscretisation(double t0, std::string hstr): TimeDiscretisation{t0,0.}
{
  mpf_set_str(_hgmp, hstr.c_str(), 10);
  mpf_set_d(_t0gmp, t0);
}

siconos::simulation::TimeDiscretisation::TimeDiscretisation(unsigned int nSteps, double t0, double T)
    : TimeDiscretisation{t0, (T - t0) / nSteps}
{}

siconos::simulation::TimeDiscretisation::TimeDiscretisation(const std::vector<double>& tk)
  : TimeDiscretisation{tk.at(0), 0.}
{
  _tkV = tk;
}


// Copy constructor
siconos::simulation::TimeDiscretisation::TimeDiscretisation(const TimeDiscretisation& td):
TimeDiscretisation{td.getTkVector()}
{
  if (td.hGmp()) {
    mpf_set(_hgmp, *td.currentTimeStep());
    _h = 0.;
  }
  else if (td.hConst()) {
    _h = td._h;
  }
  else {
    _h = 0.;
  }
}

// --- Destructor ---
siconos::simulation::TimeDiscretisation::~TimeDiscretisation()
{
  if (!_tkV.empty()) _tkV.clear();

  mpf_clear(_hgmp);
  mpf_clear(_tkp1);
  mpf_clear(_tk);
  mpf_clear(_t0gmp);
}

void siconos::simulation::TimeDiscretisation::setTkVector(const std::vector<double>& newTk)
{
  _tkV.clear();
  _tkV = newTk;
}

void siconos::simulation::TimeDiscretisation::setT0(double val)
{
  _t0 = val;
  if (_h == 0.0) mpf_set_d(_t0gmp, val);
  if (!_tkV.empty())
    THROW_EXCEPTION(
        "siconos::simulation::TimeDiscretisation::setT0 must be called only when the TimeDiscretisation is with a "
        "constant h");
}

double siconos::simulation::TimeDiscretisation::timeStep(unsigned int k)
{
  if (_tkV.empty()) {
    if (_h > 0.)
      return _h;
    else {
      mpf_mul_ui(_tkp1, _hgmp, k + 1);
      mpf_mul_ui(_tk, _hgmp, k);
      mpf_add(_tk, _tk, _t0gmp);
      mpf_add(_tkp1, _tkp1, _t0gmp);
      return mpf_get_d(_tkp1) - mpf_get_d(_tk);
    }
  }
  else
    return _tkV.at(k + 1) - _tkV.at(k);
}

double siconos::simulation::TimeDiscretisation::getTk(unsigned int indx)
{
  if (_tkV.empty()) {
    if (_h > 0.)
      return _t0 + _h * indx;
    else {
      mpf_mul_ui(_tk, _hgmp, indx);
      mpf_add(_tk, _tk, _t0gmp);
      return mpf_get_d(_tk);
    }
  }
  else
    return _tkV.at(indx);
}

// --- Other functions ---
void siconos::simulation::TimeDiscretisation::display() const
{
  std::cout << "====> Time Disretisation :\n";
  std::cout << " the current timestep is " << _h << "\n";
  std::cout << "====\n";
}
