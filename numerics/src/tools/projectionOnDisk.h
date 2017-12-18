/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#ifndef ProjectionOnDisk_H
#define ProjectionOnDisk_H

/*!\file projectionOnDisk.h
 * \brief functions to project on a cylinder
 */

#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** projectionOnDisk Projection onto the positive Disk of radius R  \f$  \{ r,  0 \sqrt(r_1^2+r_2^2) \geq R \} \f$
  \param[in,out] r the vector to be projected
  \param[in] R the radius of the cone
  */
  void projectionOnDisk(double* r, double  R);

  /** projectionOnGeneralDisk Projection onto the positive Disk of radius R  \f$  \{ r, 0 \sqrt(r_1^2+r_1^2) \geq R \} \f$
  \param[in,out] r the vector to be projected
  \param[in] R the radius of the cone
  \param[in] dim dimension of the cylinder
  */
  void projectionOnGeneralDisk(double* r, double  R, int dim);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
