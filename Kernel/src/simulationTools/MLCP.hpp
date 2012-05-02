/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
/*! \file MLCP.hpp
\brief Linear Complementarity Problem formulation and solving
*/

#ifndef MLCP_H
#define MLCP_H

#include "LinearOSNS.hpp"
#define MLCP_NB_BLOCKS 30
/** Formalization and Resolution of a Mixed Linear Complementarity Problem (MLCP)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * \section MLCPintro Aim of the MLCP class
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
 *  - a NonSmoothSolver, used to define a solver and its parameters (connected to SolverOptions structure of Numerics)
 *
 *  A MLCP is connected to a simulation that handles a NonSmoothDynamicalSystem and its Topology. \n
 *  IndexSets from simulation are used to know which constraints (Interaction) are active or not. \n
 *
 * \b Construction:
 *   - XML reading (inputs = xml node with tag "OneStepNSProblem" and a SP::Simulation)
 *   - Constructor from data (inputs = Simulations*, id, SP::NonSmoothSolver) - The solver is optional.
 * Main functions:
 *
 * \b Main functions:
 *  - formalization of the problem: computes A,B,C,D,a,b using the set of "active" Interactions from the simulation and \n
 *  the interactionBlock-matrices saved in the field interactionBlocks.\n
 *  Functions: initialize(), computeInteractionBlock(), preCompute()
 *  - solving of the problem: function compute(), used to call solvers from Numerics through \n
 * the mlcp_driver() interface of Numerics.
 *  - post-treatment of data: set values of y/lambda variables of the active Interaction (ie Interactions) using \n
 *  ouput results from the solver (u,v,z); function postCompute().
 *
 *
 */
class MLCP : public LinearOSNS
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MLCP);


  /** n is the number of equality */
  int _n;

  /** m is the size of the complementarity conditions */
  int _m;

  int _curBlock;

  /** The MLCP instance */
  MixedLinearComplementarityProblem _numerics_problem;

  /** default constructor (private)
   */
  //  MLCP(); ambiguous

public:

  /** xml constructor
  *  \param SP::OneStepNSProblemXML : the XML linked-object
  */
  MLCP(SP::OneStepNSProblemXML);

  /** constructor from data
  *  \param Solver* pointer to object that contains solver algorithm and formulation \n
  *  (optional, default = NULL => read .opt file in Numerics)
  *  \param String: id of the problem (default = "unamed")
  */
  MLCP(const int newNewNumericsSolverId = SICONOS_MLCP_ENUM);

  /** destructor
  */
  ~MLCP() {};
  virtual void computeOptions(SP::Interaction inter1, SP::Interaction inter2);
  virtual void updateInteractionBlocks();
  virtual void computeAllInteractionBlocks();
  // --- n ---
  /** get the value of n,
  *  \return int
  */
  inline int getn() const
  {
    return _n;
  }

  // --- numerics MLCP ---
  /** get the pointer on the Numerics MLCP,
  *  \return SP::MixedLinearComplementarityProblem
  */
  inline SP::MixedLinearComplementarityProblem getNumericsMLCP()
  {
    return createSPtrMixedLinearComplementarityProblem(_numerics_problem);
  }


  /** */
  virtual void reset();

  /** compute extra-diagonal interactionBlock-matrix
   *  \param an edge descriptor
   */
  virtual void computeInteractionBlock(const InteractionsGraph::EDescriptor&);

  /** compute diagonal Interaction block
   * \param a vertex descriptor
   */
  virtual void computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor&);

  virtual void preCompute(double time);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
  *  \param double : current time
  *  \return int, information about the solver convergence.
  */
  int compute(double);

  /** initialize
   */
  void initialize(SP::Simulation sim);


  /** print the data to the screen
  */
  virtual void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static MLCP* convert(OneStepNSProblem* osnsp);

};

TYPEDEF_SPTR(MLCP);
#endif // MLCP_H
