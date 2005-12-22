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
#ifndef NLGSSOLVER_H
#define NLGSSOLVER_H

#include "Solver.h"
#include "NLGSSolverXML.h"

/** \class NLGSSolver
 *  \brief to define parameters for NLGS solver
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) December 19, 2005
 *
 *  This class inherits from Solver one.
 *
 */

class NLGSSolver : public Solver
{

private:

  /** maximum iterations number */
  unsigned int maxIter;

  /** algorithm tolerance */
  double tolerance;

  /** \fn void setSolvingMethod();
   *  \brief Function to fill structure solvingMethod fields
   */
  void setSolvingMethod();

public:

  /** \fn NLGSSolver(const std::string& = DefaultSolvingForm, const unsigned int& = DefaultAlgoMaxIter, const double & = DefaultAlgoTolerance)
   *  \brief constructor with parameters values
   *  \param string : solving formalisation
   *  \param unsigned int : max iterations number
   *  \param double : tolerance
   */
  NLGSSolver(const std::string& = DefaultSolvingForm, const unsigned int& = DefaultAlgoMaxIter, const double & = DefaultAlgoTolerance);

  /** \fn NLGSSolver(const NLGSSolver&)
   *  \brief copy constructor
   *  \param a NLGSSolver
   */
  NLGSSolver(const NLGSSolver&);

  /** \fn NLGSSolver(SolverXML*)
   *  \brief constructor with XML object of the NLGSSolver
   *  \param a pointer to SolverXML
   */
  NLGSSolver(SolverXML*);

  /** \fn ~NLGSSolver()
   *  \brief destructor
   */
  ~NLGSSolver();

  // GETTERS/SETTERS

  /** \fn const unsigned int getMaxIter() const
   *  \brief get maximum iterations number
   *  \return an unsigned int
   */
  inline const unsigned int getMaxIter() const
  {
    return maxIter;
  };

  /** \fn void setMaxIter(const unsigned int&)
   *  \brief set maximum iterations number
   *  \param an unsigned int
   */
  inline void setMaxIter(const unsigned int& newVal)
  {
    maxIter = newVal;
  };

  /** \fn const double getTolerance() const
   *  \brief get tolerance
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

  /** \fn NLGSSolver* convert (Solver* solv)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Solver* : the solver that must be converted
   * \return a pointer on the NLGS Solver if it is of the right type, NULL otherwise
   */
  NLGSSolver* convert(Solver*);
};

#endif // NLGSSolver_H
