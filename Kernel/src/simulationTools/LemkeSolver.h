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
#ifndef LEMKESOLVER_H
#define LEMKESOLVER_H

#include "Solver.h"
#include "LemkeSolverXML.h"

//#include "LemkeSolverXML.h"

/** \class LemkeSolver
 *  \brief to define parameters for Lemke solver
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) December 19, 2005
 *
 *  This class inherits from Solver one.
 *
 */
class LemkeSolver : public Solver
{

private:

  /** maximum iterations number */
  unsigned int maxIter;

  /** \fn void setSolvingMethod();
   *  \brief Function to fill structure solvingMethod fields
   */
  void setSolvingMethod();

public:

  /** \fn LemkeSolver(const unsigned int & max)
   *  \brief constructor with parameters values
   *  \param string : solving formalisation type (optional -> else default value)
   *  \param unsigned int : max iterations number (optional)
   */
  LemkeSolver(const std::string& = DefaultSolvingForm, const unsigned int& = DefaultAlgoMaxIter);

  /** \fn LemkeSolver(const LemkeSolver&)
   *  \brief copy constructor
   *  \param a LemkeSolver
   */
  LemkeSolver(const LemkeSolver&);

  /** \fn LemkeSolver(SolverXML*)
   *  \brief constructor with XML object of the LemkeSolver
   *  \param a pointer to SolverXML
   */
  LemkeSolver(SolverXML*);


  /** \fn ~LemkeSolver()
   *  \brief destructor
   */
  ~LemkeSolver();

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

  /** \fn void display()
   *  \brief display solver data
   */
  void display() const;

  /** \fn LemkeSolver* convert (Solver* solv)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Solver* : the solver that must be converted
   * \return a pointer on the Lemke Solver if it is of the right type, NULL otherwise
   */
  LemkeSolver* convert(Solver*);
};

#endif // LemkeSolver_H
