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
/*! \file MixedComplementarityConditionNSL.hpp

*/
#ifndef MIXEDCOMPLEMENTARITYCONDITIONNSLAW_H
#define MIXEDCOMPLEMENTARITYCONDITIONNSLAW_H

#include "NonSmoothLaw.hpp"

namespace siconos::modeling {
/** Complementarity NonSmoothLaw
 *
 **/
class MixedComplementarityConditionNSL : public NonSmoothLaw {
 private:
  ACCEPT_SERIALIZATION(MixedComplementarityConditionNSL);

  /** default constructor
   */
  MixedComplementarityConditionNSL() = default;
  siconos::algebra::Index _equalitySize{0};

 public:
  /** basic constructor
   *
   *  \param newSize size of the non smooth law
   *  \param equalitySize size of the equality relation
   */
  MixedComplementarityConditionNSL(siconos::algebra::Index newSize,
                                   siconos::algebra::Index equalitySize)
      : NonSmoothLaw(newSize + equalitySize), _equalitySize{equalitySize} {};

  /** Destructor */
  ~MixedComplementarityConditionNSL() noexcept = default;

  /** print the data to the screen
   */
  inline void display() const override {};

  /**\return the number of equality present in the MLCP */
  inline auto equalitySize() { return _equalitySize; };

  // visitors hook
  virtual void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // MIXEDCOMPLEMENTARITYCONDITIONNSLAW_H
