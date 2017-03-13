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

#ifndef AVI_CST_H
#define AVI_CST_H
/*!\file AVI_cst.h
  \brief Constants to define the list of available AVI solvers. See the solver list \ref aviSolversList
*/
/**\enum AVI_SOLVER
   Each SICONOS_AVI_XXX refers to number of the solver XXX for an AVI. See the solver list \ref aviSolversList
 */
enum AVI_SOLVER
{
  SICONOS_AVI_CAOFERRIS = 800,
  SICONOS_AVI_PATHAVI   = 801,
};

extern const char* const   SICONOS_AVI_CAOFERRIS_STR;
extern const char* const   SICONOS_AVI_PATHAVI_STR;

#endif
