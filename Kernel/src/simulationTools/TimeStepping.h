/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#ifndef TIMESTEPPING_H
#define TIMESTEPPING_H

#include "Strategy.h"
/** \class TimeStepping
 *  \brief Specific strategy, using time stepping schemes.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.2.
 *  \date (Creation) Apr 26, 2004
 *
 */
class TimeStepping : public Strategy
{
public:

  /** \fn TimeStepping(Model * = NULL)
   *  \brief Default constructor
   *  \param a pointer to the model that owns this strategy. NULL Model leads to exception
   */
  TimeStepping(Model* = NULL);

  /** \fn TimeStepping(Model&)
   *  \brief constructor from Model => avoid this function, prefer the one with Model*
   *  \param a Model.
   */
  TimeStepping(Model&);

  /** \fn TimeStepping(StrategyXML*, Model*)
   *  \brief constructor with XML object for TimeStepping
   *  \param StrategyXML* : the XML object corresponding
   *  \param Model* : the Model which contains the Strategy
   */
  TimeStepping(StrategyXML*, Model*);

  ~TimeStepping();

  /** \fn void initialize()
   *  \brief executes the complete initialisation of Strategy (OneStepIntegrators, OneStepNSProblem, TImediscretisation) with the XML Object
   */
  void initialize();

  /** \fn void run()
   *  \brief run the simulation, from t0 to T
   */
  void run();

  /** \fn void computeOneStep()
   *  \brief run one step of the simulation
   */
  void computeOneStep();

  /** \fn TimeStepping* convert (Strategy* str)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Strategy* : the Strategy which must be converted
   * \return a pointer on the Strategy if it is of the right type, NULL otherwise
   */
  static TimeStepping* convert(Strategy* str);
};

#endif // TIMESTEPPING_H
