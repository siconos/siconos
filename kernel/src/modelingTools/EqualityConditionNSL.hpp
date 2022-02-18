/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
/*! \file EqualityConditionNSL.hpp

*/
#ifndef EQUALITYCONDITIONNSLAW_H
#define EQUALITYCONDITIONNSLAW_H

#include "NonSmoothLaw.hpp"

/** Equality NonSmoothLaw
 *
 **/
class EqualityConditionNSL : public NonSmoothLaw {
private:
  // serialization hooks
  ACCEPT_SERIALIZATION(EqualityConditionNSL);

  /** default constructor
   */
  EqualityConditionNSL() = default;

public:
  /** basic constructor
   *  \param size of the non smooth law
   */
  EqualityConditionNSL(unsigned int size);

  /** Destructor */
  ~EqualityConditionNSL();

  /** display the data of the NonSmoothLaw on the standard output
   *
   */
  void display() const override{};

  // Visitors hook
  ACCEPT_STD_VISITORS();
};

TYPEDEF_SPTR(EqualityConditionNSL)

#endif // EQUALITYCONDITIONNSLAW_H
