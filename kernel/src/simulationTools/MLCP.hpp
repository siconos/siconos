/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include <mlcp_cst.h>

#define MLCP_NB_BLOCKS 200
/** Formalization and Resolution of a Mixed Linear Complementarity Problem (MLCP)
 
   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) Apr 26, 2004
 
  \section MLCPintro Aim of the MLCP class
 
  This class is devoted to the formalization and the resolution of the
  Mixed Linear Complementarity Problem (MLCP) defined by :
    \f[
  0 =  Au + Cv + a
  \f]
  \f[
  z =  Du + Bv + b
  \f]
  \f[
  v \geq 0, z \geq 0,  z^{T} v =0
  \f]
  where
     - \f$ u \in R^{n} \f$ \f$ v \in R^{m} \f$  and \f$z \in R^{m} \f$ are the unknowns,
     - \f$ a \in R^{n} \f$ and \f$ b \in R^{m} \f$
     - \f$ A \in R^{n \times n } \f$
     - \f$ B \in R^{m \times m } \f$
     - \f$ C \in R^{n \times m } \f$
     - \f$ D \in R^{m \times n } \f$
 
   The MLCP main components are:
   - a problem (variables A,B,C,D,a,b and size of the problem), which directly corresponds to the MixedLinearComplementarityProblem structure of Numerics
   - the unknowns u,v and z
   
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

public:

  /** constructor from data
  *  \param numericsSolverId id of Numerics solver
  *  (optional, SICONOS_MLCP_ENUM the enumerative solver)
  */
  MLCP(int numericsSolverId = SICONOS_MLCP_ENUM);

  /** destructor
  */
  virtual ~MLCP() {reset();};

  /** compute equalities/inequalities sizes and set corresponding values in
      numerics problem
      \param inter1 Interaction used to get a non-smooth law and the constraints sizes.
      \param inter2 another interaction, not used indeed (?)
  */
  virtual void computeOptions(SP::Interaction inter1, SP::Interaction inter2);

  /** Update blocks used to compute M matrix.
   */
  virtual void updateInteractionBlocks();

  /** get the number of equality constraints,
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

  /** Reninitialize numerics driver.
   */
  virtual void reset();

  /** compute extra-diagonal interactionBlock-matrix
   *  \param ed an edge descriptor
   */
  virtual void computeInteractionBlock(const InteractionsGraph::EDescriptor& ed);

  /** compute diagonal Interaction block
   * \param vd a vertex descriptor
   */
  virtual void computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd);

  /** Pre compute 
   * \param time current time
   * \return bool
   */
  virtual bool preCompute(double time);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
  *  \param time current time
  *  \return int, information about the solver convergence.
  */
  int compute(double time);

  /** initialize
   * \param sim the Simulation
   */
  void initialize(SP::Simulation sim);


  /** print the data to the screen
  */
  virtual void display() const;

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

#endif // MLCP_H
