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
#ifndef LCP_AS_CONVEXQP_H
#define LCP_AS_CONVEXQP_H

/*!\file LinearComplementarityProblem_as_ConvexQP.h
  \brief Definition of a structure to handle with LCP problems.
*/

#include "NumericsFwd.h"
#include "SiconosConfig.h"

/** \struct LinearComplementarityProblem_as_ConvexQP LinearComplementarityProblem_as_ConvexQP.h
 *
 */
struct LinearComplementarityProblem_as_ConvexQP
{
  /* the ConvexQP associated with the FC3D problem */
  ConvexQP * cqp;
  /* the FC3D associated with the ConvexQP  */
  LinearComplementarityProblem * lcp;
  /* the SolverOptions that might be used to pass some numerical parameters */
  SolverOptions * options;
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void Projection_ConvexQP_LCP(void *cqpIn, double *x, double *PX);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
