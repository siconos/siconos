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
#include "LinearComplementaritySystemsNSDS.hpp"
#include "Interaction.hpp"
#include "Relation.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "ComplementarityConditionNSL.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "debug.h"


using namespace RELATION;

// --- CONSTRUCTORS/DESTRUCTOR ---

//  constructor
LinearComplementaritySystemsNSDS::LinearComplementaritySystemsNSDS(double t0, double T, SP::SiconosVector x0,
                                                                   SP::SimpleMatrix A,  SP::SimpleMatrix B,
                                                                   SP::SimpleMatrix C,  SP::SimpleMatrix D,
                                                                   SP::SiconosVector a,  SP::SiconosVector b): NonSmoothDynamicalSystem(t0,T)
{
  _ds.reset(new FirstOrderLinearTIDS(x0, A));
  if (a)
  {
    _ds-> setbPtr(a);
  }
  insertDynamicalSystem(_ds);

  _relation.reset(new FirstOrderLinearTIR(C, B));

  // todo: check sizes
  if (D)
  {
    _relation->setDPtr(D);
  }
  if (b)
  {
    _relation->setePtr(b);
  }
  _nslaw.reset(new ComplementarityConditionNSL(C->size(0)));
  _interaction.reset(new Interaction(_nslaw, _relation));

  link(_interaction, _ds);

  display();
};
