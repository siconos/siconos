/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#ifndef NCP_CST_H
#define NCP_CST_H

/*!\file NCP_cst.h
  \brief Constants to define the list of available NCP solvers. See the solver list \ref ncpSolversList
*/

/**\enum NCP_SOLVER
   Each SICONOS_NCP_XXX refers to number of the solver XXX for NCP. See the solver list \ref ncpSolversList
 */
enum NCP_SOLVER
{
  SICONOS_NCP_NEWTON_FBLSA = 900,
  SICONOS_NCP_NEWTON_MINFBLSA = 901,
  SICONOS_NCP_PATHSEARCH = 902,
  SICONOS_NCP_PATH = 903
};


extern const char* const  SICONOS_NCP_NEWTON_FBLSA_STR;
extern const char* const  SICONOS_NCP_NEWTON_MINFBLSA_STR;
extern const char* const  SICONOS_NCP_PATHSEARCH_STR;
extern const char* const  SICONOS_NCP_PATH_STR;

#endif
