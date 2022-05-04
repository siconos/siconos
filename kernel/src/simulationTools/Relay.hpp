/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
/*! \file Relay.hpp
  \brief Linear Complementarity Problem formulation and solving
*/

#ifndef Relay_H
#define Relay_H

#include "LinearOSNS.hpp"

#include <relay_cst.h>
#include <RelayProblem.h>
TYPEDEF_SPTR(RelayProblem)

/**
   Formalization and Resolution of a Linear Complementarity Problem (Relay)

   This class is devoted to the formalization and the resolution of the
   Relay NonSmooth problems.

   \f[
   w =  q + M z
   \f]
   \f[
   w \geq 0, z \geq 0,  z^{T} w =0
   \f]
   where
   - \f$ w \in R^{n} \f$  and \f$ z \in R^{n} \f$ are the unknowns,
   - \f$ M \in R^{n \times n } \f$  and \f$ q \in R^{n} \f$
   
   \todo : add "recover" function to start from old values of z and w.
   \todo : review this introduction ...
*/
class Relay : public LinearOSNS
{

protected:
  
  ACCEPT_SERIALIZATION(Relay);


  /** contains the vector lb (lower bounds) of a Relay system */
  SP::SiconosVector _lb;

  /** contains the vector ub (upper bounds) of a Relay system */
  SP::SiconosVector _ub;

  /** contains the numerics proble for Relay system */
  SP::RelayProblem _numerics_problem;

  /** nslaw effects : visitors experimentation
   */

  struct _BoundsNSLEffect;
  friend struct _BoundsNSLEffect;



public:

  /** constructor from numerics solver id
   *
   *  \param numericsSolverId id of numerics solver, default =  SICONOS_RELAY_AVI_CAOFERRIS
   */
  Relay(int numericsSolverId = SICONOS_RELAY_AVI_CAOFERRIS);

  /** constructor from a pre-defined solver options set
   *
   *  \param options the options set
   */
  Relay(SP::SolverOptions options);

  /** destructor
   */
  ~Relay(){};

  // --- lb ---
  /** get the value of lb, the   lower bounds of the Relay system
   *
   *  \return the vector of lower bounds
   */
  inline const SiconosVector& getLb() const
  {
    return *_lb;
  }

  /** get lb, the lower bounds of the Relay system
   *
   *  \return the vector of lower bounds
   */
  inline SP::SiconosVector lb() const
  {
    return _lb;
  }

  /** set lb to pointer newPtr
   *
   *  \param newLb new lower bound
   */
  inline void setLb(SP::SiconosVector newLb)
  {
    _lb = newLb;
  }


  // --- ub ---
  /** get the value of ub, the  upper bounds of the Relay system
   *
   *  \return the vector of upper bounds
   */
  inline const SiconosVector& getUb() const
  {
    return *_ub;
  }

  /** get lb, the lower bounds of the Relay system
   *
   *  \return the vector of upper bounds
   */
  inline SP::SiconosVector ub() const
  {
    return _ub;
  }

  /** set ub to pointer newPtr
   *
   *  \param newUb new upper bound
   */
  inline void setUb(SP::SiconosVector newUb)
  {
    _ub = newUb;
  }

  void initialize(SP::Simulation sim);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
   *
   *  \param time current time
   *  \return information about the solver convergence.
   */
  int compute(double time);

  /** Check the compatibility fol the nslaw with the targeted OSNSP */
  bool checkCompatibleNSLaw(NonSmoothLaw& nslaw);


  /** print the data to the screen
   */
  void display() const;

};

#endif // Relay_H
