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

#ifndef SA_SCAL_HPP
#define SA_SCAL_HPP

/** multiplication of a matrix by a scalar, B = a*A (init = true) or B += a*A (init = false)
 *  \param a a double
 *  \param A a SiconosMatrix
 *  \param[in,out] B a SiconosMatrix
 *  \param init a bool
 */
void scal(double a, const SiconosMatrix& A, SiconosMatrix& B, bool = true);



#endif
