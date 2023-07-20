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

#include "NM_MA57.h"
#ifdef WITH_MA57

#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"


/*#define DEBUG_MESSAGES*/
#include "siconos_debug.h"

void NM_MA57_free(void* p)
{
  NSM_linear_solver_params* params = (NSM_linear_solver_params*) p;
  if (params->linear_solver_data){
    LBL_Data * lbl =  (LBL_Data *)params->linear_solver_data;
    FILE * logfile = lbl->ma57->logfile;
    LBL_Finalize(lbl);
    fclose(logfile);
  }
  //free(params->linear_solver_data);
  params->linear_solver_data = NULL;
}

#endif
