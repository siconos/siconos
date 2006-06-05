/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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

/** \class OneStepNSProblem
 *  \brief It's the part of the Strategy which solve the Interactions
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.2.0.
 *  \date (Creation) Apr 26, 2004
 *
 *
 * This is an abstract class, that provides an interface to define a non smooth problem:
 *   -> a formulation (ie the way the problem is written)
 *   -> a solver (algorithm and solving formulation, that can be different from problem formulation)
 *   -> routines to compute the problem solution.
 *
 * The available problem formulation, given by derived classes, are:
 *  - LCP
 *  - FrictionContact2D and 3D
 *  - QP
 *  - Relay
 *
 *  See Solver class or Numerics documentation for details on algorithm name and parameters.
 *
 */

#ifndef ONESTEPNSPROBLEM_H
#define ONESTEPNSPROBLEM_H

#include "Strategy.h"
#include "Interaction.h"
#include "EqualityConstraint.h"
#include "SiconosMatrix.h"
#include "SimpleVector.h"
#include "OneStepNSProblemXML.h"
#include "Topology.h"
#include "SiconosConst.h"
#include "SiconosNumerics.h"
#include "check.h"
#include "Solver.h"
#include <iostream>
#include <vector>
#include <string>

class Strategy;
class Interaction;
class EqualityConstraint;
class SiconosMatrix;
class OneStepNSProblemXML;
class Solver;
class OneStepNSProblem
{

protected:

  /** type of the OneStepNSProblem (LCP ...) */
  std::string nspbType;

  /** size of the problem to solve */
  unsigned int dim;

  /** all the Interactions known by the OneStepNSProblem */
  InteractionsSet OSNSInteractions;

  /** all the Equality Constraints known by the OneStepNSProblem */
  std::vector< EqualityConstraint* > ecVector;

  /** map that links each interaction with the corresponding diagonal block*/
  std::map< Interaction* , SiconosMatrix*>  diagonalBlocksMap  ;

  /** map that links each interaction with the corresponding extra-diagonal blocks
      map < InteractionA * , map < InteractionB* , blockMatrixAB > >
      Interaction A and B are coupled through blockMatrixAB.
  */
  std::map< Interaction* , std::map<Interaction *, SiconosMatrix*> >  extraDiagonalBlocksMap  ;

  /** map that links each interaction with a list of indexes, giving the indexes
      corresponding to the effective relations AND their derivatives */
  std::map< Interaction* , std::vector<unsigned int> > blockIndexesMap ;

  /** Solver for Non Smooth Problem*/
  Solver* solver;

  /** bool to check whether solver has been allocated inside the class or not */
  bool isSolverAllocatedIn;

  /** link to the strategy that owns the NSPb */
  Strategy *strategy;

  /** the XML object linked to the OneStepNSProblem to read XML data */
  OneStepNSProblemXML* onestepnspbxml;

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** \fn OneStepNSProblem()
   *  \brief default constructor
   */
  OneStepNSProblem();

  /** \fn OneStepNSProblem(OneStepNSProblemXML*, Strategy*=NULL)
   *  \brief xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Strategy *: the strategy that owns the problem (optional)
   */
  OneStepNSProblem(OneStepNSProblemXML*, Strategy * = NULL);

  /** \fn OneStepNSProblem(Strategy*, Solver*)
   *  \brief constructor from data
   *  \param Strategy *: the strategy that owns this problem
   *  \param Solver *: pointer on object that contains solver algorithm definition (optional)
   */
  OneStepNSProblem(Strategy * , Solver* = NULL);

  /** \fn OneStepNSProblem()
   *  \brief destructor
   */
  virtual ~OneStepNSProblem();

  // --- GETTERS/SETTERS ---

  /** \fn inline const string getType() const
   *  \brief to get the type of the OneStepNSProblem
   *  \return string
   */
  inline std::string getType() const
  {
    return nspbType;
  }

  /** \fn inline void setType(const string&)
   *  \brief set the type of the OneStepNSProblem
   *  \param: string
   */
  inline void setType(const std::string & newVal)
  {
    nspbType = newVal;
  }

  /** \fn const int getN() const
   *  \brief get dimension of the problem
   *  \return an unsigned ing
   */
  inline const unsigned int getDim() const
  {
    return dim;
  }

  /** \fn void setDim(const int&)
   *  \brief set the value of dim
   *  \param an unsigned int
   */
  inline void setDim(const unsigned int& newVal)
  {
    dim = newVal;
  }

  /** \fn const InteractionsSet getInteractions()
    *  \brief get the set of Interactions associated with the NS Problem
    *  \return an InteractionsSet
    */
  inline const InteractionsSet getInteractions() const
  {
    return OSNSInteractions;
  };

  /** \fn void setInteractionst(const InteractionsSet&)
   *  \brief set the Interaction list of this NS Problem
   *  \param an InteractionsSet
   */
  void setInteractions(const InteractionsSet&);

  /** \fn Interaction* getInteractionPtr(const int&)
   *  \brief get a specific Interaction
   *  \param int the position of a specific Interaction in the vector of Interaction
   *  \return a pointer on Interaction
   */
  Interaction* getInteractionPtr(const unsigned int&);

  /** \fn vector< EqualityConstraint* > getEqualityConstraints()
   *  \brief get the the vector of EqualityConstraint
   *  \return a vector stl
   */
  inline const std::vector< EqualityConstraint* > getEqualityConstraints() const
  {
    return ecVector;
  }

  /** \fn void setEqualityConstraints(vector< EqualityConstraint* >)
    *  \brief set the vector of EqualityConstraint
    *  \param vector<EqualityConstraint*> : a vector stl
    */
  inline void setEqualityConstraints(const std::vector< EqualityConstraint* >& newVec)
  {
    ecVector = newVec;
  }

  /** \fn Strategy* getStrategyPtr()
   *  \brief get the Strategy
   *  \return a pointer on Strategy
   */
  inline Strategy* getStrategyPtr() const
  {
    return strategy;
  }

  /** \fn void setStrategyPtr(Strategy*)
   *  \brief set the Strategy of the OneStepNSProblem
   *  \param: a pointer on Strategy
   */
  inline void setStrategy(Strategy* str)
  {
    strategy = str;
  }

  /** \fn Solver* getSolverPtr()
   *  \brief get the Solver
   *  \return a pointer on Solver
   */
  inline Solver* getSolverPtr() const
  {
    return solver;
  }

  /** \fn void setSolverPtr(Solver*)
   *  \brief set the Solver of the OneStepNSProblem
   *  \param: a pointer on Solver
   */
  void setSolverPtr(Solver*);

  /** \fn inline OneStepNSProblemXML* getOneStepNSProblemXML()
   *  \brief get the OneStepNSProblemXML
   *  \return a pointer on OneStepNSProblemXML
   */
  inline OneStepNSProblemXML* getOneStepNSProblemXML() const
  {
    return onestepnspbxml;
  }

  /** \fn inline void setOneStepNSProblemXML(OneStepNSProblemXML* osnspb)
   *  \brief set the OneStepNSProblemXML
   *  \param a pointer on OneStepNSProblemXML
   */
  inline void setOneStepNSProblemXML(OneStepNSProblemXML* osnspb)
  {
    onestepnspbxml = osnspb;
  }

  // --- OTHER FUNCTIONS ---

  /** \fn void addInteraction(Interaction*)
   *  \brief add an Interaction to the OneStepNSProblem
   *  \param Interaction* : the Interaction to add to the vector of Interaction
   */
  void addInteraction(Interaction*);

  /** \fn void initialize()
   *  \brief initialize the problem(compute topology ...)
   */
  virtual void initialize() = 0;

  /** \fn void computeEffectiveOutput();
   *  \brief compute variables indexMax, effectiveOutput and effectiveSizeOutput of the topology of the nsds
   */
  void computeEffectiveOutput();

  /** \fn void nextStep(void)
   *  \brief prepares the problem for the next time step
   *  \exception to be defined
   */
  void nextStep();

  /** \fn void updateInput()
   *  \brief compute r thanks to lambda
   */
  void updateInput();

  /** \fn void updateOutput(void)
   *  \brief compute output for all the interactions
   */
  void updateOutput();

  /** \fn void compute(const double&)
   *  \brief make the computation so solve the NS problem
   *  param double : current time
   */
  virtual void compute(const double&) = 0;

  /** \fn void saveNSProblemToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveNSProblemToXML() = 0;

  /** \fn void isOneStepNsProblemComplete()
   *  \brief check is all data required by non smooth problem are present
   *  \return : a bool
   */
  bool isOneStepNsProblemComplete() const;

  /** \fn void check_solver(const int& info) const
   *  \brief return exception and message if solver failed
   *  \param: output from solve_... (Numerics routine)
   */
  void check_solver(const int&) const;

};

#endif // ONESTEPNSPROBLEM_H
