/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
/*! \file MLCPProjectOnConstraints.hpp
\brief Linear Complementarity Problem formulation and solving
*/

#ifndef MLCPProjectOnConstraints_H
#define MLCPProjectOnConstraints_H

#include "MLCP.hpp"
#include "OSNSMatrixProjectOnConstraints.hpp"
/** Formalization and Resolution of a Mixed Linear Complementarity Problem (MLCP)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.3.0.
 *  \date (Creation) 01, 2011
 *
 * This class is devoted to the formalization and the resolution of the
 * Mixed Linear Complementarity Problem (MLCP) for the specific problem
 * of the projection onto the constraints in Mechanics
 *
 */
class MLCPProjectOnConstraints : public MLCP
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MLCPProjectOnConstraints);
  
  /** ?? */
  double _alpha;

  /** disabled or enabled projection On Equality (or Unilateral) for unilateral constraints */
  bool _doProjOnEquality;

  bool  _useMassNormalization;

public:

  /** compute the number of inequality and equality for a given tuple of Interactions
      update the global number of equality(_n) and inequality (_m)
      set up _numerics_problem parameters (blocksRows and blocksIsComp )
      \param inter1 first interaction considered
      \param inter2 second interaction
  */
  virtual void computeOptions(SP::Interaction inter1, SP::Interaction inter2);


  /** constructor from data
   \param numericsSolverId solver id
   \param alpha alpha parameter value
   */
  MLCPProjectOnConstraints(int numericsSolverId = SICONOS_MLCP_ENUM, double alpha = 1.0);

  /** destructor
  */
  ~MLCPProjectOnConstraints() {};

  /** 
      \return alpha value
  */
  double alpha()
  {
    return _alpha;
  };

  /** setter for alpha
      \param[in] newval new value for alpha parameter
  */
  void setAlpha(double newval)
  {
    _alpha = newval;
  };

  inline void setDoProjOnEquality(bool v)
  {
    _doProjOnEquality = v;
  }



  /** Display the set of blocks for  a given indexSet
      \param indexSet the graph of interactions
   */
  void displayBlocks(SP::InteractionsGraph indexSet);

  /** print the data to the screen
  */
  void display() const;

  virtual void initOSNSMatrix();

  /** compute interactionBlocks if necessary (this depends on the type of
   * OSNS, on the indexSets ...)
   */
  virtual void updateInteractionBlocks();

  /** compute interactionBlocks if necessary (this depends on the type of
   * OSNS, on the indexSets ...)
   */
  virtual void updateInteractionBlocksOLD();

  /** compute diagonal Interaction block
      \param vd a vertex (interaction) descriptor
  */
  virtual void computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd);

  /** compute diagonal Interaction block
      \param vd a vertex (interaction) descriptor
  */
  virtual void computeInteractionBlock(const InteractionsGraph::EDescriptor& vd);

  /** To compute a part of the "q" vector of the OSNS
      \param vd vertex (interaction) which corresponds to the considered block
      \param pos the position of the first element of yOut to be set
  */
  virtual void computeqBlock(InteractionsGraph::VDescriptor& vd, unsigned int pos);

  /** compute vector q
   *  \param time the current time
   */
  void computeq(double time);
  
  /** post-treatment for  MLCPProjectOnConstraints
   */
  virtual void postCompute();

  /** post-treatment for  MLCPProjectOnConstraints for LagrangianR
      \param inter the considered interaction
      \param pos interaction position in the global vector
  */
  virtual void postComputeLagrangianR(SP::Interaction inter, unsigned int pos);

  /** post-treatment for  MLCPProjectOnConstraints for NewtonEulerR
      \param inter the considered interaction
      \param pos interaction position in the global vector
   */
  virtual void postComputeNewtonEulerR(SP::Interaction inter, unsigned int pos);

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

#endif // MLCPProjectOnConstraints_H
