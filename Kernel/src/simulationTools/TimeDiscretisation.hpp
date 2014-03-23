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
/*! \file TimeDiscretisation.hpp

  \brief A time-discretisation scheme for a given interval.
 */
#ifndef TIMEDISCRETISATION_H
#define TIMEDISCRETISATION_H

#include "TimeDiscretisationXML.hpp"
#include <vector>
#include <gmp.h>

/** A time discretisation scheme

    \author SICONOS Development Team - copyright INRIA
    \version 3.0.0.
    \date (Creation) Apr 26, 2004

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

  /** the XML object linked to the TimeDiscretisation to read XML data */
  SP::TimeDiscretisationXML _timeDiscretisationXML;

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

  /** constructor with XML
      \param tdXML the XML object corresponding
      \param t0 initial time
      \param T final time
  */
  TimeDiscretisation(SP::TimeDiscretisationXML tdXML, double t0, double T);

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

  /** get the TimeDiscretisationXML of the TimeDiscretisation
   *  \return a pointer on the TimeDiscretisationXML of the TimeDiscretisation
   */
  inline SP::TimeDiscretisationXML timeDiscretisationXML() const
  {
    return _timeDiscretisationXML;
  }

  /** set the TimeDiscretisationXML of the TimeDiscretisation
   *  \param timediscrxml a SP::TimeDiscretisationXML
   */
  inline void setTimeDiscretisationXMLPtr(SP::TimeDiscretisationXML timediscrxml)
  {
    _timeDiscretisationXML = timediscrxml;
  }

  // --- OTHER FUNCTIONS ---
  /** print the discretisation data to the screen
   */
  void display() const;

  // --- XML Functions ---

  /** saves the TimeDiscretisation to the XML tree
   *
   */
  void saveTimeDiscretisationToXML();
};

//DEFINE_SPTR(TimeDiscretisation)

#endif // TIMEDISCRETISATION_H
