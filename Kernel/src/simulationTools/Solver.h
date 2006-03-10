/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
 *  \version 1.1.3.
 *  \date (Creation) December 19, 2005
 *
 *  It provides an interface to define a solver for a given solving formalisation (the way the problem is written to be solved, ie the non
 *  smooth problem type). The solver characteristics are:
 *    - the way the problem is written, defined by the class type name of the OneStepNSProblem that owns this solver
 *    - the name of the solver, ie the numerical algorithm that is used, among those provide by Numerics package
 *    - a structure method from Numerics (see related doc for more details)
 *    - some specific parameters, which are required or not, depending on the used algorithm.
 *
 *  Available solvers algorithms are (and their required parameters):
 *   - Lemke (maxIter) => warning: not fully implemented in Numerics, do not use!
 *   - LexicoLemke (maxIter)
 *   - NLGS (Non Linear Gauss Seidel) (tolerance, maxIter)
 *   - NSQP (Non Smooth Quadratic Programming) (tolerance)
 *   - QP   (tolerance)
 *   - CPG  (Conjugate Projected Gradient) (tolerance, maxIter)
 *   - Latin (tolerance, maxIter, searchDirection)
 *
 *  Notes:
 *          - parameter normType is never used, but present for future developments.
 *          - it is Numerics work to test whether non-smooth problem accept a solver or not, nothing is done in that sense in the Kernel.
 *
 */

#include "SiconosConst.h"
#include "SiconosNumerics.h"
#include "RuntimeException.h"
#include<string>
#include<iostream>

// Default values for non smooth problem solver
const std::string DEFAULT_SOLVER = "NLGS";
const double DEFAULT_TOL = 0.0001;
const unsigned int DEFAULT_ITER = 101;
const std::string DEFAULT_NORMTYPE = "max";
const double DEFAULT_SEARCHDIR = 0.6;

class SolverXML;

class Solver
{
protected:

  /** non smooth problem formulation
   *  This corresponds to the type of the OneStepNSProblem
   *  that owns this solver */
  std::string nonSmoothPbType;

  /** name of the used solver algorithm
      This corresponds to the derived class name*/
  std::string solverAlgorithmName;

  /** Numerics object that contains information and structure of various
      solving methods. The present class is dedicated to the filling of this
      object. */
  method* solvingMethod;

  // Warning: some of the following data may be useless, depending on the solver type.
  // See full documentation for more details.

  /** maximum iterations number */
  unsigned int maxIter;

  /** algorithm tolerance */
  double tolerance;

  /** */
  std::string normType;

  /**  */
  double searchDirection;

  /** \fn virtual void setSolvingMethod() = 0;
   *  \brief Function to fill structure solvingMethod fields
   */
  void setSolvingMethod();

  /** \fn Solver(const Solver&)
   *  \brief default constructor
   */
  Solver();

public:

  /** \fn Solver(const std::string& = DefaultAlgoName, const std::string& = DefaultSolvingForm)
   *  \brief constructor with the value of the Solver attributes
   *  \param string: non smooth problem type (LCP ...)
   *  \param string: solver name (optional)
   *  \param unsigned int   : MaxIter (optional) required if a solver is given
   *  \param double : Tolerance (optional) -> for NLGS, Gcp, Latin
   *  \param string : NormType (optional) -> for NLGS, Gcp, Latin
   *  \param double : SearchDirection (optional) -> for Latin
   */
  Solver(const std::string&, const std::string& = DEFAULT_SOLVER, const unsigned int & = DEFAULT_ITER,
         const double & = DEFAULT_TOL, const std::string & = DEFAULT_NORMTYPE, const double & = DEFAULT_SEARCHDIR);

  /** \fn Solver(const Solver&)
   *  \brief copy constructor
   *  \param a Solver to be copied
   */
  Solver(const Solver&);

  /** \fn Solver(SolverXML*)
   *  \brief constructor with XML object of the Solver
   *  \param a pointer to SolverXML
   *  \param string: non smooth problem type
   */
  Solver(SolverXML*, const std::string&);

  /** \fn ~Solver()
   *  \brief destructor
   */
  ~Solver();

  // GETTERS/SETTERS

  /** \fn const string getNonSmoothPbType() const
   *  \brief get the solver algorithm name
   *  \return a string
   */
  inline const std::string getNonSmoothPbType() const
  {
    return nonSmoothPbType;
  };

  /** \fn void setNonSmoothPbType(const string&)
   *  \brief set the solver algorithm name
   *  \param a string
   */
  inline void setNonSmoothPbType(const std::string& newVal)
  {
    nonSmoothPbType = newVal;
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

  /** \fn const double getSearchDirection() const
   *  \brief get searchDirection
   *  \return an double
   */
  inline const double getSearchDirection() const
  {
    return searchDirection;
  };

  /** \fn void setSearchDirection(const double&)
   *  \brief set searchDirection
   *  \param an double
   */
  inline void setSearchDirection(const double& newVal)
  {
    searchDirection = newVal;
  };

  /** \fn const std::string getNormType() const
   *  \brief get normType
   *  \return a string
   */
  inline const std::string getNormType() const
  {
    return normType;
  };

  /** \fn void setNormType(const std::string&)
   *  \brief set normType
   *  \param an std::string
   */
  inline void setNormType(const std::string& newVal)
  {
    normType = newVal;
  };

  /** \fn void saveSolverToXML()
   *  \brief copy the data of the Solver into the XML tree
   */
  void saveSolverToXML();

  /** \fn void display()
   *  \brief display solver data
   */
  void display() const ;
};

#endif // Solver_H
