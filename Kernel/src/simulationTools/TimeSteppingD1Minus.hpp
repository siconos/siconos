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

/*!\file
 * TimeSteppingD1Minus simulation
 */

#ifndef TIMESTEPPINGD1MINUS_H
#define TIMESTEPPINGD1MINUS_H

#include "Simulation.hpp"

/** TimeSteppingD1Minus Timestepping Strategy
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) September 16, 2011
 *
 *  see Schindler/Acary : Timestepping Schemes for Nonsmooth Dynamics Based
 *  on Discontinuous Galerkin Methods: Definition and Outlook
 */

class TimeSteppingD1Minus : public Simulation
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(TimeSteppingD1Minus);

  /** default constructor */
  TimeSteppingD1Minus() {}

protected:
  /** initialisation specific to TimeSteppingD1Minus for OneStepNSProblem */
  virtual void initOSNS();

public:

  /** constructor with the time-discretisation
   *  \param pointer to a TimeDiscretisation
   *  \param number of non smooth problem
   */
  TimeSteppingD1Minus(SP::TimeDiscretisation, int nb);

  /** destructor */
  ~TimeSteppingD1Minus();

  /* type name */
  virtual std::string typeName()
  {
    return Type::name(*this);
  }

  /** updateIndexSet using current y and lambda values of interactions
   *  \param unsigned int: number of the set to be updated
   *  0 : ALL interactions (NEVER)
   *  1 : ACTIVE interactions
   *  2 : STAYING ACTIVE interactions
   */
  virtual void updateIndexSet(unsigned int);

  /** update input, state and output of DynamicalSystems
   *  \param level to be updated for input
   */
  virtual void update(unsigned int);

  /** run the simulation, from t0 to T */
  virtual void run();

  /** step from current event to next event of EventsManager */
  virtual void advanceToEvent();

  /** compute residu */
  void computeResidu();

  /** integrate DynamicalSystems taking not into account non-smooth part */
  void computeFreeState();

  /** encapsulates an operation of dynamic casting
   *  needed by Python interface
   *  \param simulation which must be converted
   *  \return pointer to the Simulation if it is of the right type, NULL otherwise
   */
  static TimeSteppingD1Minus* convert(Simulation* str);

  /** visitors hook */
  ACCEPT_STD_VISITORS();
};

DEFINE_SPTR(TimeSteppingD1Minus);

#endif // TIMESTEPPINGD1MINUS_H
