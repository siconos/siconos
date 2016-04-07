/* Siconos-Kernel, Copyright INRIA 2005-2015
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
/*! \file AVI.hpp
  \brief Affine Variational Inequalities formulation
*/

#ifndef AVI_H
#define AVI_H

#include "LinearOSNS.hpp"

#include <AVI_cst.h>
#include <AffineVariationalInequalities.h>
TYPEDEF_SPTR(AffineVariationalInequalities)

/** Formalization and Resolution of a Linear Complementarity Problem (AVI)
 
   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) Apr 26, 2004
 
  \section AVIintro Aim of the AVI class
 
  This class is devoted to the formalization and the resolution of
  Affine variational Inequalities (AVI): given a polytopic set \f$P\f$, \f$M\in R^{p\times p}\f$ and \f$q\in R^p\f$,
  \f[
  \text{find z}\in P\text{such that for all x}\in P\quad \langle x-z, Mz+q\rangle \geq 0.
  \f]
  \todo : add "recover" function to start from old values of z and w.
*/
class AVI : public LinearOSNS
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(AVI);

  /** contains the numerics problem for the AVI system */
  SP::AffineVariationalInequalities _numerics_problem;

public:

  /** constructor from data
   *  \param numericsSolverId id of numerics solver
   */
  AVI(int numericsSolverId = SICONOS_AVI_CAOFERRIS);

  /** destructor
   */
  ~AVI();

  void initialize(SP::Simulation sim);
  virtual void setSolverId(int solverId);

  /** Compute the unknown z and update the Interaction (y and lambda)
   *  \param time current time
   *  \return information about the solver convergence.
   */
  int compute(double time);

  /** print the data to the screen
   */
  void display() const;

};

#endif // AVI_H
