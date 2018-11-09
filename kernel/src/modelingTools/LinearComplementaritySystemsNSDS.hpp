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
/*! \file NonSmoothDynamicalSystem.hpp
 * \brief container for DynamicalSystem and Interaction
 */
#ifndef LinearComplementaritySystemsNSDS_H
#define LinearComplementaritySystemsNSDS_H

#include "SiconosPointers.hpp"
#include "Topology.hpp"
#include "DynamicalSystem.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "ComplementarityConditionNSL.hpp"
#include "SiconosFwd.hpp"
#include "FirstOrderLinearTIR.hpp"


/** the LinearComplementaritySystemsNSDS_H inherits frim NDSD
    for a direct instanciation of a LCS
*/
class LinearComplementaritySystemsNSDS: public NonSmoothDynamicalSystem
{


private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LinearComplementaritySystemsNSDS);


  /* a first order linear TI dynamical systems */
  SP::FirstOrderLinearTIDS _ds;
  /* a first order linear TI relation */
  SP::FirstOrderLinearTIR _relation;
  /* a complementarity condition */
  SP::ComplementarityConditionNSL _nslaw;
  /* an interaction*/
  SP::Interaction _interaction;


public:

  /** default constructor
   */
  LinearComplementaritySystemsNSDS();

  /** constructor with t0 and T
   * \param t0 initial time
   * \param T final time
   */
  LinearComplementaritySystemsNSDS(double t0, double T,  SP::SiconosVector x0,
                                   SP::SimpleMatrix A,  SP::SimpleMatrix B,
                                   SP::SimpleMatrix C,  SP::SimpleMatrix D,
                                   SP::SiconosVector a,  SP::SiconosVector b);

  /** destructor
   */
  ~LinearComplementaritySystemsNSDS();

  // --- GETTERS/SETTERS ---

  SP::FirstOrderLinearTIDS ds()
  {
    return _ds;
  };

  SP::FirstOrderLinearTIR relation()
  {
    return _relation;
  };

  SP::ComplementarityConditionNSL nslaw()
  {
    return _nslaw;
  };
  SP::Interaction interaction()
  {
    return _interaction;
  };


};


#endif
