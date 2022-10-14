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

siconos::simulation::TimeDiscretisation::TimeDiscretisation(double t0, double h)
    : default_time_step_(h), t0_(t0)
{
  mpf_init(_hgmp);
  mpf_init(_tkp1);
  mpf_init(_tk);
  mpf_init(_t0gmp);
  // Set timeStep method
  timeStep = [this](std::size_t = 0) { return default_time_step_; };
  getTk = [this](std::size_t indx) { return t0_ + default_time_step_ * indx; };
}

siconos::simulation::TimeDiscretisation::TimeDiscretisation(double t0, std::string hstr)
    : TimeDiscretisation{t0, 0.}
{
  gmp_is_on_ = true;
  mpf_set_str(_hgmp, hstr.c_str(), 10);
  mpf_set_d(_t0gmp, t0);
  timeStep = [this](std::size_t k = 0) {
    mpf_mul_ui(_tkp1, _hgmp, k + 1);
    mpf_mul_ui(_tk, _hgmp, k);
    mpf_add(_tk, _tk, _t0gmp);
    mpf_add(_tkp1, _tkp1, _t0gmp);
    return mpf_get_d(_tkp1) - mpf_get_d(_tk);
  }; // FP : why don't we only return _hgmp supposed to be constant ???

  getTk = [this](std::size_t indx) {
    mpf_mul_ui(_tk, _hgmp, indx);
    mpf_add(_tk, _tk, _t0gmp);
    return mpf_get_d(_tk);
  };
}

siconos::simulation::TimeDiscretisation::TimeDiscretisation(unsigned int nSteps, double t0,
                                                            double T)
    : TimeDiscretisation{t0, (T - t0) / nSteps}
{
}

siconos::simulation::TimeDiscretisation::TimeDiscretisation(const std::vector<double>& tk)
    : time_instants_{tk}
{
  step_is_constant_ = false;
  assert(tk.size() > 1 && "Please provide a vector which size is at least 2.");
  t0_ = tk[0];
  timeStep = [this](std::size_t k) { return time_instants_.at(k + 1) - time_instants_.at(k); };
  // 'at' to throw exception if the requested index does not exist.

  getTk = [this](std::size_t indx) { return time_instants_.at(indx); };
}

// Copy constructor
siconos::simulation::TimeDiscretisation::TimeDiscretisation(const TimeDiscretisation& td)
    : TimeDiscretisation{0, 0.}
{
  if (gmp_is_on_) {
    mpf_set(_hgmp, *td.currentTimeStep());
    default_time_step_ = 0.;
  }
  else if (step_is_constant_) {
    default_time_step_ = td.default_time_step_;
  }
  else {
    default_time_step_ = 0.;
  }
  t0_ = td.getT0();
  time_instants_ = td.time_instants_;
  timeStep = td.timeStep;
  getTk = td.getTk;
}

// --- Destructor ---
siconos::simulation::TimeDiscretisation::~TimeDiscretisation() noexcept
{
  mpf_clear(_hgmp);
  mpf_clear(_tkp1);
  mpf_clear(_tk);
  mpf_clear(_t0gmp);
}

// --- Other functions ---
void siconos::simulation::TimeDiscretisation::display() const
{
  if (hConst()) {
    std::cout << "====> Fixed time-step time discretisation :" << std::endl;
    std::cout << " the current timestep is " << timeStep(0) << "\n";
  }
  else
    std::cout << "====> Variable time-step time discretisation.\n";
  std::cout << "====\n";
}
