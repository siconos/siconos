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
/*! \file TimeDiscretisation.hpp

  \brief A time-discretisation scheme for a given interval.
 */
#ifndef TIMEDISCRETISATION_H
#define TIMEDISCRETISATION_H

#include <vector>
#include <gmp.h>
#include "SiconosFwd.hpp"
#include "SiconosSerialization.hpp" // for ACCEPT_SERIALIZATION

typedef std::vector<double> TkVector;

/** A time discretisation scheme

    A TimeDiscretisation object is used to discretized a given time interval. \n
    TimeDiscretisation are used:
    - in the simulation, as a user-input to discretized [t0,T]
    - in Sensor and Actuator, to define the set of time instants where the sensor or actuator
    must operate.

    A TimeDiscretisation is defined with a starting time (t0), a time step size (h, non necessarily constant),
    the number of the current time step (k). \n
    The time instant values are saved in a vector tk. Depending on the way of construction of the TimeDiscretisation,
    all or only current and next times are saved in tk. The difference is just a question of saving memory. \n

    Note that the TimeDiscretisation is not linked to the Model. It's up to the user to check that the way he builds his time-discretisation fits with the t0 and T given in the model.

    \section tdMfunc Main functions:
    - setCurrentTimeStep(), to set current h. This value will be used for all future time steps, until next change.
    - increment(), shift to next time step (increment k, and shift t[k] and t[k+1])
    - currentTime(), return t[k]

    \section tdBuild Construction

    - input = the complete vector tk. This defines t0, T, number of time steps and time step size
    (which is not necessarily constant). In this case, the whole vector is saved in the memory.
    - inputs = number of time steps, t0 and T.
    size h is computed with t0,T and nSteps. Only two values are saved: t[k] and t[k+1] = t[k] + h.
    h can be changed at any time.
    - inputs = h and t0. Only two values are saved: t[k] and t[k+1] = t[k] + h.
    h can be changed at any time.
    - inputs = t0 and h.  Only two values are saved: t[k] and t[k+1] = t[k] + h.
    h can be changed at any time.

 */
class TimeDiscretisation
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(TimeDiscretisation);


  /** Default value for the time step (tk+1 - tk) */
  double _h;

  /** vector of time values at each step (=> size = n1+n2+1 - Default size = 2 - Max size= nSteps+1) */
  TkVector _tkV;

  /** Origin of time*/
  double _t0;

  /** Timestep stored as mpf_t, for high precision computations */
  mpf_t _hgmp;

  /** Time at t_{k+1}, in mpf_t, used to compute a good timestep */
  mpf_t _tkp1;

  /** Time at t_{k+1}, in mpf_t, used to compute a good timestep */
  mpf_t _tk;

   /** Starting time, used to compute a good timestep */
  mpf_t _t0gmp;

  /** default constructor (private => no copy nor pass-by value)
   */
  TimeDiscretisation(); 

  /** Assignment Operator (private => forbidden)
   * \param td unused
   * \return TimeDiscretisation&
   */
  TimeDiscretisation& operator =(const TimeDiscretisation& td);

  /** Get the origin of time t0
   * \return the origin of time
   */
  inline double getT0() const
  {
    return _t0;
  }

public:

  // --- Straightforward constructors ---

  /** constructor with tk vector of instant times values.
   *  \param newTk a TkVector describing the discretisation
   */
  TimeDiscretisation(const TkVector& newTk);

  /** constructor with the size of the default time step and t0
   *  \param t0 initial time value
   *  \param h the time step
   */
  TimeDiscretisation(double t0, double h);

  /** constructor with the number of steps, t0 and T
   *  \param nSteps the number of steps
   *  \param t0 initial time value
   *  \param T the final time
   */
  TimeDiscretisation(unsigned int nSteps, double t0, double T);

  /** Constructor with the size of the default timestep and t0.
   * It creates a TimeDiscretisation using GMP for all its computation
   *  \param t0 initial time value
   *  \param str the time step in form of a string
   */
  TimeDiscretisation(double t0, const std::string& str);

  /** Copy constructor
   * \param td the TimeDiscretisation to copy
   */
  TimeDiscretisation(const TimeDiscretisation& td);

  // Destructor
  ~TimeDiscretisation();

  // --- GETTERS/SETTERS ---

  /** get the timestep \f$t_{k+1} - t_k\f$
   * \param k the index of the timestep
   * \return the current time step
   */
  double currentTimeStep(const unsigned int k);

  /** get the timestep in gmp format
   * \return a pointer to the timestep in mpf_t format
   */
  inline const mpf_t* currentTimeStep() const
  {
    return &_hgmp;
  };

  /** determine whether the timestep is constant
   * \return true if the timestep is constant
   */
  inline bool hConst() const
  {
    return _tkV.empty() ? true : false;
  };

  /** determine whether the TimeDiscretisation is using GMP
   * \return true if the TimeDiscretisation is using GMP
   */
  inline bool hGmp() const
  {
    return ((_h == 0.0) && (_tkV.empty()))  ? true : false;
  }

  /** get the value of tk at step k
   * \param indx the step
   * \return a double
   */
  double getTk(const unsigned int indx);

  /** get the TkVector _tkV
   *  \return a reference to the TkVector _tkV
   */
  inline const TkVector& getTkVector() const { return _tkV; };

  /** set the TkVector _tkV
   *  \param newTk the new value for _tkV
   */
  void setTkVector(const TkVector& newTk);

  /** change t0 before the simulation starts (useful for delays)
   *  \param val the new value for t0
   */
  void setT0(double val);

  // --- OTHER FUNCTIONS ---
  /** print the discretisation data to the screen
   */
  void display() const;

};

#endif // TIMEDISCRETISATION_H
