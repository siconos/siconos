/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
/*! \file Equality.hpp
  \brief Linear Complementarity Problem formulation and solving
*/

#ifndef Equality_H
#define Equality_H

#include "LinearOSNS.hpp"

/** Formalization and Resolution of a Linear Complementarity Problem (Equality)
 
  \section Equalityintro Aim of the Equality class
 
  This class is devoted to the formalization and the resolution of the
  Linear system (Equality) defined by :
  \f$
  0 = w =  q + M z
  \f$
  where

     - \f$ w \in R^{n} \f$  and \f$z \in R^{n} \f$ are the unknowns,
     - \f$ M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 
   The Equality main components are:
   - a problem (variables M,q and size of the problem), which directly corresponds to the LinearComplementarityProblem structure of Numerics
   - the unknowns z and w
 
*/
class Equality : public LinearOSNS
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Equality);

public:

  /** constructor from data
   *  \param numericsSolverId id of numerics solver, default = 0
   */
  Equality(int numericsSolverId = 0);

  /**  constructor from a pre-defined solver options set.
       \param options, the options set, 
       \rst
       see :ref:`problems_and_solvers` for details.
       \endrst
  */
  Equality(SP::SolverOptions options);
  
  /** destructor
   */
  ~Equality() {};

  /** initialize
   * \param sim the simulation
   */
  void initialize(SP::Simulation sim);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
   *  \param time double : current time
   *  \return int information about the solver convergence.
   */
  int compute(double time);

  /** Build or reinit M and the NumericsProblem*/
  virtual void updateM();

   /* Check the compatibility fol the nslaw with the targeted OSNSP */
  bool checkCompatibleNSLaw(NonSmoothLaw& nslaw);

  /** print the data to the screen
   */
  void display() const;
};

#endif // Equality_H
