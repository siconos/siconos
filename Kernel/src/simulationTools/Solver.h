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
#ifndef SOLVER_H
#define SOLVER_H

#include "SolverXML.h"

/** \class Solver
 *  \brief to define a solver for a Non smooth problem.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) December 19, 2005
 *
 *  !!! This class is an abstract one !!!
 *
 *  It provides an interface to define a solving formalisation (the way the problem is written to be solved),
 *  a solver algorithm (which corresponds to one of the derived class) and its
 *  required parameters (fields of the derived class).
 *  The method structure from Numerics is embedded in the present class.
 *
 * Various solving formalisations are:
 *   - lcp ( "LcpSolving")
 *   -
 *
 * Available solvers algorithms are:
 *   - Lemke
 *   - LexicoLemke
 *   - NLGS (Non Linear Gauss Seidel)
 *   - NSQP (Non Smooth Quadratic Programming
 *   - QP
 *   - CPG  (Conjugate Projected Gradient)
 *   - Latin
 *
 */

#include "SiconosConst.h"
#include "SiconosNumerics.h"
#include "RuntimeException.h"
#include<string>
#include<iostream>
class SolverXML;

class Solver
{
protected:

  /** name of the solvingFormalisation to use
      Warning: this variable may be different from NSProblem type */
  std::string solvingFormalisation ;

  /** name of the used solver algorithm
      This corresponds to the derived class name*/
  std::string solverAlgorithmName;

  /** Numerics object that contains information and structure of various
      solving methods. The present class is dedicated to the filling of this
      object. */
  method* solvingMethod;

  /** \fn virtual void setSolvingMethod() = 0;
   *  \brief Function to fill structure solvingMethod fields
   */
  virtual void setSolvingMethod() = 0;

public:

  /** \fn Solver(const std::string& = DefaultAlgoName, const std::string& = DefaultSolvingForm)
   *  \brief constructor with the value of the Solver attributes
   *  \param string: solver algorithm name (optional)
   *  \param string: solving formalisation type (optional)
   *  Note that without any parameters, this is the default constructor
   */
  Solver(const std::string& = DefaultAlgoName, const std::string& = DefaultSolvingForm);

  /** \fn Solver(const Solver&)
   *  \brief copy constructor
   *  \param a Solver to be copied
   */
  Solver(const Solver&);

  /** \fn Solver(SolverXML*)
   *  \brief constructor with XML object of the Solver
   *  \param a pointer to SolverXML
   */
  Solver(SolverXML*);

  /** \fn ~Solver()
   *  \brief destructor
   */
  virtual ~Solver();

  // GETTERS/SETTERS

  /** \fn const string getSolvingFormalisation() const
   *  \brief get the solving formalisation type
   *  \return a string
   */
  inline const std::string getSolvingFormalisation() const
  {
    return solvingFormalisation;
  };

  /** \fn void setSolvingFormalisation(const string&)
   *  \brief set the solving formalisation type
   *  \param a string
   */
  inline void setSolvingFormalisation(const std::string& newVal)
  {
    solvingFormalisation = newVal;
  };

  /** \fn const string getSolverAlgorithmName() const
   *  \brief get the solver algorithm name
   *  \return a string
   */
  inline const std::string getSolverAlgorithmName() const
  {
    return solverAlgorithmName;
  };

  /** \fn void setSolverAlgorithmName(const string&)
   *  \brief set the solver algorithm name
   *  \param a string
   */
  inline void setSolverAlgorithmName(const std::string& newVal)
  {
    solverAlgorithmName = newVal;
  };

  /** \fn method* getSolvingMethodPtr() const
   *  \brief get method (Numerics) structure
   *  \return a pointer to method
   */
  inline method* getSolvingMethodPtr() const
  {
    return solvingMethod;
  };

  /** \fn void saveSolverToXML()
   *  \brief copy the data of the Solver into the XML tree
   */
  virtual void saveSolverToXML();

  /** \fn void display()
   *  \brief display solver data
   */
  virtual void display() const = 0 ;

  /** \fn Solver* convert (Solver* solv)
    *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
    *  \param Solver* : the solver that must be converted
    * \return a pointer on the XXXSolver if it is of the right type, XXX being the solver name
    * NULL otherwise
    */
  virtual Solver* convert(Solver*) = 0;
};

#endif // Solver_H
