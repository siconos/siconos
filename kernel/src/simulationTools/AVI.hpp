/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
/*! \file AVI.hpp
  \brief Affine Variational Inequalities formulation
*/

#ifndef AVI_H
#define AVI_H

#include "LinearOSNS.hpp"

#include <AVI_cst.h>
#include <AffineVariationalInequalities.h>
TYPEDEF_SPTR(AffineVariationalInequalities)

/** Formalization and Resolution of an Affine Variational Inequality (AVI)
 
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
