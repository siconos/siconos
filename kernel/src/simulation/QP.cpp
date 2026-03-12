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

#include "QP.hpp"

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

int siconos::nonsmooth_formulations::QP::compute(double)
{
  THROW_EXCEPTION("siconos::nonsmooth_formulations::QP::compute not yet implemented");
  return 1;
}

void siconos::nonsmooth_formulations::QP::display() const
{
  THROW_EXCEPTION("siconos::nonsmooth_formulations::QP::compute not yet implemented");
}

void siconos::nonsmooth_formulations::QP::setQ(const siconos::algebra::SiconosMatrix& newValue)
{
  *_Q = newValue;
}
void siconos::nonsmooth_formulations::QP::setP(const siconos::algebra::SiconosVector& newValue)
{
  *_p = newValue;
}
