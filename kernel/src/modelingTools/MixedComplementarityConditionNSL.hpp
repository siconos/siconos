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
/*! \file MixedComplementarityConditionNSL.hpp

*/
#ifndef MIXEDCOMPLEMENTARITYCONDITIONNSLAW_H
#define MIXEDCOMPLEMENTARITYCONDITIONNSLAW_H

#include "NonSmoothLaw.hpp"

#include "SiconosPointers.hpp"


/** Complementarity NonSmoothLaw
 *
 **/
class MixedComplementarityConditionNSL : public NonSmoothLaw
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MixedComplementarityConditionNSL);

  /** default constructor
   */
  MixedComplementarityConditionNSL() {};
  unsigned int _equalitySize;

public:
  /** basic constructor
   *  \param newSize size of the non smooth law
   *  \param equalitySize size of the equality relation
   */
  MixedComplementarityConditionNSL(unsigned int newSize, unsigned int equalitySize);

  /** Destructor */
  ~MixedComplementarityConditionNSL();


  /** print the data to the screen
  */
  inline void display()const {};

  /** get the number of equality present in the MLCP
   *  \return an unsigned int
   */
  inline unsigned int equalitySize()
  {
    return _equalitySize;
  };

  /** Visitors hook
   */
  ACCEPT_STD_VISITORS();
};

#endif // MIXEDCOMPLEMENTARITYCONDITIONNSLAW_H
