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

#include "AlgebraTools.hpp"
#include "SiconosMatrix.hpp"
#include <boost/numeric/ublas/io.hpp>
#include "expm.hpp"

namespace Siconos
{
namespace algebra
{
namespace tools
{

void expm(SiconosMatrix& A, SiconosMatrix& Exp, bool computeAndAdd)
{
  // Implemented only for dense matrices.
  // Note FP : Maybe it works for others but it has not been
  // tested here --> to be done
  // Do not work with sparse.
  A.resetFactorizationFlags();
  Exp.resetFactorizationFlags();
  assert(Exp.num() == Siconos::DENSE || A.num() == Siconos::DENSE);
  if(computeAndAdd)
    *Exp.dense() += expm_pad(*A.dense());
  else
    *Exp.dense() = expm_pad(*A.dense());
}
} // namespace tools
} // namespace algebra
} // namespace Siconos
