/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#ifndef ProjectionOnCylinder_H
#define ProjectionOnCylinder_H

#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** projectionOnCylinder Projection onto the positive Cylinder of radius R  \f$  \{ r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq R \} \f$
  \param[in,out] r the vector to be projected
  \param[in] R the radius of the cone
  */
  void projectionOnCylinder(double* r, double  R);
  
  /** projectionOnGeneralCylinder Projection onto the positive Cylinder of radius R  \f$  \{ r, r_1 \geq 0, 0 \sqrt(r_2^2+r_3^2) \geq R \} \f$
  \param[in,out] r the vector to be projected
  \param[in] R the radius of the cone
  \param[in] dim dimension of the cylinder
  */
  void projectionOnGeneralCylinder(double* r, double  R, int dim);
  
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
