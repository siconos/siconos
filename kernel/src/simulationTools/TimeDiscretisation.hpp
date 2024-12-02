/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
/*! \file TimeDiscretisation.hpp

  A time-discretisation scheme for a given interval.
 */
#ifndef TIMEDISCRETISATION_H
#define TIMEDISCRETISATION_H

#include <gmp.h>

#include <limits>
#include <vector>
#include <functional>

#include "SiconosSerialization.hpp"  // for ACCEPT_SERIALIZATION

namespace siconos::simulation {

/**
   A time discretisation scheme

   A TimeDiscretisation object is used to discretized a given time interval.
   TimeDiscretisation are used:
   - in the simulation, as a user-input to discretized [t0,T]
   - in Sensor and Actuator, to define the set of time instants where the sensor or actuator
   must operate.

   A TimeDiscretisation is defined with a starting time (t0), a time step size (h, non
   necessarily constant), the number of the current time step (k). The time instant values are
   saved in a vector tk. Depending on the way of construction of the TimeDiscretisation, all or
   only current and next times are saved in tk. The difference is just a question of saving
   memory.

   Note that the TimeDiscretisation is not linked to the Model. It's up to the user to check
   that the way he builds his time-discretisation fits with the t0 and T given in the model.

   Main functions:
   - setCurrentTimeStep(), to set current h. This value will be used for all future time steps,
   until next change.
   - increment(), shift to next time step (increment k, and shift t[k] and t[k+1])
   - currentTime(), return t[k]

   Construction

   - input = the complete vector tk. This defines t0, T, number of time steps and time step
   size (which is not necessarily constant). In this case, the whole vector is saved in the
   memory.
   - inputs = number of time steps, t0 and T.
   size h is computed with t0,T and nSteps. Only two values are saved: t[k] and t[k+1] = t[k] +
   h. h can be changed at any time.
   - inputs = h and t0. Only two values are saved: t[k] and t[k+1] = t[k] + h.
   h can be changed at any time.
   - inputs = t0 and h.  Only two values are saved: t[k] and t[k+1] = t[k] + h.
   h can be changed at any time.

*/
class TimeDiscretisation {
 private:
  ACCEPT_SERIALIZATION(TimeDiscretisation);

  // Function proto used to get current time step or current time value
  using TimeFunction = std::function<double(std::size_t)>;

  /** Default value for the time step (tk+1 - tk) */
  double default_time_step_{0.};

  /** vector of time values at each step (=> size = n1+n2+1 - Default size = 2 - Max size=
   * nSteps+1) */
  std::vector<double> time_instants_ = {};

  /** Origin of time*/
  double t0_{std::numeric_limits<double>::quiet_NaN()};

  /** Check if time-step is constant or not */
  bool step_is_constant_{true};

  // Note FP : gmp stuff should be moved into a derived class ? Or CRTP ?
  // Or removed since never used ...
  /** Check if gmp is on **/
  bool gmp_is_on_{false};

  /** Timestep stored as mpf_t, for high precision computations */
  mpf_t _hgmp;

  /** Time at t_{k+1}, in mpf_t, used to compute a good timestep */
  mpf_t _tkp1;

  /** Time at t_{k+1}, in mpf_t, used to compute a good timestep */
  mpf_t _tk;

  /** Starting time, used to compute a good timestep */
  mpf_t _t0gmp;

  // Rule of five
  TimeDiscretisation() = delete;
  TimeDiscretisation(TimeDiscretisation&&) = delete;
  TimeDiscretisation& operator=(const TimeDiscretisation&) = delete;
  TimeDiscretisation& operator=(TimeDiscretisation&&) = delete;

  /** Get the origin of time t0
   * \return the origin of time
   */
  inline double getT0() const { return t0_; }

 public:
  /** constructor with the vector of instant times values. To be used only for variable
   * time-step time discretisations.
   *
   *  \param newTk a vector of double
   */
  TimeDiscretisation(const std::vector<double>& newTk);

  /** constructor with the size of the default time step and t0. Fixed time-step only.
   *
   *  \param t0 initial time value
   *  \param h the time step
   */
  TimeDiscretisation(double t0, double h);

  /** Constructor to build a 'gmp' time-discretisation.
   *  It creates a TimeDiscretisation using GMP for all its computation
   *
   *  \param t0 initial time value
   *  \param hstr the time step described as "MeN". M is the mantissa and N is the exponent
   */
  TimeDiscretisation(double t0, std::string hstr);

  /** constructor with the number of steps, t0 and T. Fixed time-step only.
   *
   *  \param nSteps the number of steps
   *  \param t0 initial time value
   *  \param T the final time
   */
  TimeDiscretisation(unsigned int nSteps, double t0, double T);

  /** Copy constructor
   *
   *  \param td the TimeDiscretisation to copy
   */
  TimeDiscretisation(const TimeDiscretisation& td);

  // Destructor
  ~TimeDiscretisation() noexcept;

  /** \return the timestep \f$ t_{k+1} - t_k \f$
   *
   *  \param k (default=0 for fixed time step) the index of the timestep.
   */
  TimeFunction timeStep{nullptr};

  /** get the value of tk at step kx
   *
   *  \param indx the step
   *  \return a double
   */
  TimeFunction getTk{nullptr};

  /** get the timestep in gmp format
   *
   *  \return a pointer to the timestep in mpf_t format
   */
  inline const mpf_t* currentTimeStep() const { return &_hgmp; };

  /** determine whether the timestep is constant
   *
   *  \return true if the timestep is constant
   */
  inline bool hConst() const {return step_is_constant_;};

  /** determine whether the TimeDiscretisation is using GMP
   *
   *  \return true if the TimeDiscretisation is using GMP
   */
  inline bool hGmp() const { return gmp_is_on_; }

  /** print the discretisation data to the screen
   */
  void display() const;
};
}  // namespace siconos::simulation

#endif  // TIMEDISCRETISATION_H
