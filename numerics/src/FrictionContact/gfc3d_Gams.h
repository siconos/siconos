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

/*!\file gfc3d_Gams.h
  \brief Various structure to work with GAMS
*/

#ifndef gfc3d_Gams_h
#define gfc3d_Gams_h

#include "NumericsMatrix.h"

/** Structure needed to solve the global friction contact problem */
typedef struct {
  NumericsMatrix* Mlu; /**< M in LU factored form*/
} GFC3D_Gams;


#endif
