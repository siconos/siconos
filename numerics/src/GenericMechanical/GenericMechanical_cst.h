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

#ifndef GenericMechanical_cst_h
#define GenericMechanical_cst_h

/*!\file GenericMechanical_cst.h
  \brief Constants to define the list of available GenericMechanical solvers. 
  See the solver list \ref genericSolversList
*/
/**\enum GENERIC_MECHANICAL_SOLVER
   \brief Each SICONOS_GENERIC_MECHANICAL_XXX refers to number of the solver XXX 
   for GENERIC_MECHANICAL. See the solver list \ref genericSolversList
 */
enum GENERIC_MECHANICAL_SOLVER
{
  SICONOS_GENERIC_MECHANICAL_NSGS = 2000
};

/** \enum iparam indices for generic mechanical solvers*/
enum GENERIC_MECHANICAL_IPARAM
  {
   SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED = 2,
   SICONOS_GENERIC_MECHANICAL_IPARAM_WITH_LINESEARCH = 19,
  };

enum GENERIC_MECHANICAL_DPARAM
  {
   SICONOS_DPARAM_GMP_ERROR_LS = 18,
   SICONOS_DPARAM_GMP_COEFF_LS = 19,
  };

/**\enum Possible values for iparam[GENERIC_MECHANICAL_IPARAM_ISREDUCED]  */
enum GENERIC_MECHANICAL_ISREDUCED
  {
   SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS = 0, // GS on all blocks
   SICONOS_GENERIC_MECHANICAL_SUBS_EQUALITIES = 1, // The equalities are substituated
   SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES = 2, // Equalities are assemblated in one block
   SICONOS_GENERIC_MECHANICAL_MLCP_LIKE = 3, // Try to solve like a MLCP (==> No FC3d)
  };

extern const char* const  SICONOS_GENERIC_MECHANICAL_NSGS_STR;

#endif
