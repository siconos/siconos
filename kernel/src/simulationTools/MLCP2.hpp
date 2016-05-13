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
/*! \file MLCP.hpp
\brief Linear Complementarity Problem formulation and solving
*/

#ifndef MLCP2_H
#define MLCP2_H

#include "MLCP.hpp"

/** Formalization and Resolution of a Mixed Linear Complementarity Problem (MLCP)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * \section MLCP2intro Aim of the MLCP2 class
 *
 * This class is devoted to the formalization and the resolution of the
 * Mixed Linear Complementarity Problem (MLCP) defined by :
 *  * \f[
 * 0 =  Au + Cv + a
 * \f]
 * \f[
 * z =  Du + Bv + b
 * \f]
 * \f[
 * v \geq 0, z \geq 0,  z^{T} v =0
 * \f]
 * where
 *    - \f$ u \in R^{n} \f$ \f$ v \in R^{m} \f$  and \f$z \in R^{m} \f$ are the unknowns,
 *    - \f$ a \in R^{n} \f$ and \f$ b \in R^{m} \f$
 *    - \f$ A \in R^{n \times n } \f$
 *    - \f$ B \in R^{m \times m } \f$
 *    - \f$ C \in R^{n \times m } \f$
 *    - \f$ D \in R^{m \times n } \f$
 *
 *  The MLCP main components are:
 *  - a problem (variables A,B,C,D,a,b and size of the problem), which directly corresponds to the MixedLinearComplementarityProblem structure of Numerics
 *  - the unknowns u,v and z
 *
 *  A MLCP is connected to a simulation that handles a NonSmoothDynamicalSystem and its Topology. \n
 *  IndexSets from simulation are used to know which constraints (Interaction) are active or not. \n
 *
 * \b Construction:
 *   - Constructor from data (inputs = Simulations*, id, NumericsSolverName) - The solver is optional.
 * Main functions:
 *
 * \b Main functions:
 *  - formalization of the problem: computes A,B,C,D,a,b using the set of "active" Interactions from the simulation and \n
 *  the block-matrices saved in the field blocks.\n
 *  Functions: initialize(), computeBlock(), preCompute()
 *  - solving of the problem: function compute(), used to call solvers from Numerics through \n
 * the mlcp_driver() interface of Numerics.
 *  - post-treatment of data: set values of y/lambda variables of the active Interaction using \n
 *  ouput results from the solver (u,v,z); function postCompute().
 *
 *
 */
class MLCP2 : public MLCP
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MLCP2);

  bool mFirstCall;

  /** default constuctor
   */
  MLCP2() {};

public:

  /** constructor from data
   * \param numericsSolverId solver to use (see numerics documentation for valid numbers)
  */
  MLCP2(int numericsSolverId = SICONOS_MLCP_PGS);;

  /** destructor
  */
  virtual ~MLCP2() {};


  /** To initialize the MLCP problem(computes topology ...)
     \param the simulation, owner of this OSNSPB
   */
  virtual void initialize(SP::Simulation);

  /** compute vector q
  *  \param double : current time
  */
  void computeq(double);

  /** pre-treatment for MLCP
  *  \param time the current time
  *  \return true if the computation has to be carry on, false otherwise
  */
  virtual bool preCompute(double time);

  virtual void computeInteractionBlock(SP::Interaction, SP::Interaction);
  virtual void computeDSBlock(SP::DynamicalSystem);
  virtual void computeInteractionDSBlock(SP::Interaction , SP::DynamicalSystem);
  virtual void computeDSInteractionBlock(SP::DynamicalSystem, SP::Interaction);
  virtual void updateM();
  /** Compute the unknown z and w and update the Interaction (y and lambda )
  *  \param double : current time
  *  \return int, information about the solver convergence.
  */
  int compute(double);



  /** post-treatment for MLCP
  */
  void postCompute() ;

};

#endif // MLCP2_H
