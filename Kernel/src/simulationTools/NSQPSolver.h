/* Siconos version 1.0, Copyright INRIA 2005.
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
#ifndef NSQPSOLVER_H
#define NSQPSOLVER_H

#include "Solver.h"
#include "NSQPSolverXML.h"

/** \class NSQPSolver
 *  \brief to define parameters for NSQP solver
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) December 19, 2005
 *
 *  This class inherits from Solver one.
 *
 */

class NSQPSolver : public Solver
{

private:

  /** algorithm tolerance */
  double tolerance;

  /** \fn void setSolvingMethod();
   *  \brief Function to fill structure solvingMethod fields
   */
  void setSolvingMethod();

public:

  /** \fn NSQPSolver(const std::string& = DefaultSolvingForm, const double & = DefaultAlgoTolerance)
   *  \brief constructor with parameters values
   *  \param string : solving formalisation
   *  \param double : tolerance
   */
  NSQPSolver(const std::string& = DefaultSolvingForm, const double & = DefaultAlgoTolerance);

  /** \fn NSQPSolver(const NSQPSolver&)
   *  \brief copy constructor
   *  \param a NSQPSolver
   */
  NSQPSolver(const NSQPSolver&);

  /** \fn NSQPSolver(SolverXML*)
   *  \brief constructor with XML object of the NSQPSolver
   *  \param a pointer to SolverXML
   */
  NSQPSolver(SolverXML*);

  /** \fn ~NSQPSolver()
   *  \brief destructor
   */
  ~NSQPSolver();

  // GETTERS/SETTERS

  /** \fn const double getTolerance() const
   *  \brief get tolerance algorithm
   *  \return an double
   */
  inline const double getTolerance() const
  {
    return tolerance;
  };

  /** \fn void setTolerance(const double&)
   *  \brief set tolerance algorithm
   *  \param an double
   */
  inline void setTolerance(const double& newVal)
  {
    tolerance = newVal;
  };

  /** \fn void display()
   *  \brief display solver data
   */
  void display() const;

  /** \fn NSQPSolver* convert (Solver* solv)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Solver* : the solver that must be converted
   * \return a pointer on the NSQP Solver if it is of the right type, NULL otherwise
   */
  NSQPSolver* convert(Solver*);
};

#endif // NSQPSolver_H
