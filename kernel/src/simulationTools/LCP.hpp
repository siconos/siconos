/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
/*! \file LCP.hpp
  Linear Complementarity Problem formulation and solving
*/

#ifndef LCP_H
#define LCP_H

#include "LinearOSNS.hpp"

#include <LinearComplementarityProblem.h>
#include <lcp_cst.h>
TYPEDEF_SPTR(LinearComplementarityProblem)

/** 
    Formalization and Resolution of a Linear Complementarity Problem (LCP)
    
    This class is devoted to the formalization and the resolution of the
    Linear Complementarity Problem (LCP) defined by :
    
    \f[ 

    w =  q + M z \\
    w \geq 0, z \geq 0,  z^{T} w =0
    \f]

    where
    - \f$ w \in R^{n} \f$  and \f$ z \in R^{n} \f$ are the unknowns,
    - \f$ M \in R^{n \times n } \f$  and \f$ q \in R^{n} \f$
    
    The LCP main components are:
    - a problem (variables M,q and size of the problem), which directly corresponds to the LinearComplementarityProblem structure of Numerics
    - the unknowns z and w
 
    For details regarding the available options, see Nonsmooth problems formulations and available solvers in users' guide.
*/
class LCP : public LinearOSNS
{
protected:

  ACCEPT_SERIALIZATION(LCP);

  /** Structure (for Numerics component) that describes the problem to solve */
  SP::LinearComplementarityProblem _numerics_problem;

public:

  /** constructor from numerics solver id
   *
   *  \param numericsSolverId id of numerics solver, default =  SICONOS_LCP_LEMKE
   */
  LCP(int numericsSolverId = SICONOS_LCP_LEMKE);

  /** constructor from a pre-defined solver options set
   *
   *  \param options the options set
   */
  LCP(SP::SolverOptions options);

  /** destructor */
  ~LCP() noexcept = default;

  /**
     Call numerics solver.
     
     \return int information about the solver convergence
     (output from numerics driver, linearComplementarity_driver, check numerics doc. for details).
   */
  int solve();

  /**
     Compute the unknowns z and w and update the corresponding Interactions (y and lambda )
     
     \param time : current time
     \return int, information about the solver convergence
     (output from numerics driver, linearComplementarity_driver, check numerics doc. for details).
  */
  int compute(double time) override;


   /* Check the compatibility fol the nslaw with the targeted OSNSP */
  bool checkCompatibleNSLaw(NonSmoothLaw& nslaw) override;

  ACCEPT_STD_VISITORS();


};

#endif // LCP_H
