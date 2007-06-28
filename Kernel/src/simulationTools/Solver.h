/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
/*! \file
*/
#ifndef SOLVER_H
#define SOLVER_H

#include "SolverXML.h"
#include "SiconosNumerics.h"

// Default values for non smooth problem solver
const std::string DEFAULT_SOLVER = "NLGS";
const double DEFAULT_TOL = 0.0001;
const unsigned int DEFAULT_ITER = 101;
const std::string DEFAULT_NORMTYPE = "max";
const double DEFAULT_SEARCHDIR = 0.6;
const unsigned int DEFAULT_VERBOSE = 0;
const double DEFAULT_RHO = 1.0;

class SolverXML;

/** Non Smooth Solver
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
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

  /** Verbose mode chat > 0*/
  unsigned int verbose;

  /** */
  std::string normType;

  /**  */
  double searchDirection;

  double Rho;

  /** default constructor
  */
  Solver();

public:

  /** constructor with the value of the Solver attributes
  *  \param string: non smooth problem type (LCP ...)
  *  \param string: solver name (optional)
  *  \param unsigned int   : MaxIter (optional) required if a solver is given
  *  \param double : Tolerance (optional) -> for NLGS, Gcp, Latin, RPGS
  *  \param string : NormType (optional) -> for NLGS, Gcp, Latin
  *  \param double : SearchDirection (optional) -> for Latin
  *  \param double : Rho (optional) -> for RPGS
  */
  Solver(const std::string&, const std::string& = DEFAULT_SOLVER, const unsigned int & = DEFAULT_ITER,
         const double & = DEFAULT_TOL, const unsigned int & = DEFAULT_VERBOSE, const std::string & = DEFAULT_NORMTYPE, const double & = DEFAULT_SEARCHDIR, const double & = DEFAULT_RHO);

  /** copy constructor
  *  \param a Solver to be copied
  */
  Solver(const Solver&);

  /** constructor with XML object of the Solver
  *  \param a pointer to SolverXML
  *  \param string: non smooth problem type
  */
  Solver(SolverXML*, const std::string&);

  /** destructor
  */
  ~Solver();

  // GETTERS/SETTERS

  /** get the solver algorithm name
  *  \return a string
  */
  inline const std::string getNonSmoothPbType() const
  {
    return nonSmoothPbType;
  };

  /** set the solver algorithm name
  *  \param a string
  */
  inline void setNonSmoothPbType(const std::string& newVal)
  {
    nonSmoothPbType = newVal;
  };

  /** get the solver algorithm name
  *  \return a string
  */
  inline const std::string getSolverAlgorithmName() const
  {
    return solverAlgorithmName;
  };

  /** set the solver algorithm name
  *  \param a string
  */
  inline void setSolverAlgorithmName(const std::string& newVal)
  {
    solverAlgorithmName = newVal;
  };

  /** get method (Numerics) structure
  *  \return a pointer to method
  */
  inline method* getSolvingMethodPtr() const
  {
    return solvingMethod;
  };

  /** get maximum iterations number
  *  \return an unsigned int
  */
  inline const unsigned int getMaxIter() const
  {
    return maxIter;
  };

  /** set maximum iterations number
  *  \param an unsigned int
  */
  inline void setMaxIter(const unsigned int& newVal)
  {
    maxIter = newVal;
  };

  /** get the verbose mode variable
  *  \return an unsigned int
  */
  inline const unsigned int getVerbose() const
  {
    return verbose;
  };

  /** set the verbose mode variable
  *  \param an unsigned int
  */
  inline void setVerbose(const unsigned int& newVal)
  {
    verbose = newVal;
  };

  /** get tolerance
  *  \return an double
  */
  inline const double getTolerance() const
  {
    return tolerance;
  };

  /** set tolerance algorithm
  *  \param an double
  */
  inline void setTolerance(const double& newVal)
  {
    tolerance = newVal;
  };

  /** get searchDirection
  *  \return an double
  */
  inline const double getSearchDirection() const
  {
    return searchDirection;
  };

  /** set searchDirection
  *  \param an double
  */
  inline void setSearchDirection(const double& newVal)
  {
    searchDirection = newVal;
  };

  /** get normType
  *  \return a string
  */
  inline const std::string getNormType() const
  {
    return normType;
  };

  /** set normType
  *  \param an std::string
  */
  inline void setNormType(const std::string& newVal)
  {
    normType = newVal;
  };

  /** get rho
  *  \return a double
  */
  inline const double getRho() const
  {
    return Rho;
  };

  /** set rho
  *  \param a double
  */
  inline void setRho(const double& newVal)
  {
    Rho = newVal;
  };

  /** Function to fill structure solvingMethod fields
  */
  void setSolvingMethod();

  /** copy the data of the Solver into the XML tree
  */
  void saveSolverToXML();

  /** display solver data
  */
  void display() const ;
};

#endif // Solver_H
