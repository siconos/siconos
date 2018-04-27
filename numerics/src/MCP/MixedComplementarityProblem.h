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
#ifndef MCP_PROBLEM_H
#define MCP_PROBLEM_H

/*!\file MixedComplementarityProblem.h
*/

#include "SiconosConfig.h"
#include "NumericsFwd.h"

/** type for user defined function used to compute Fmcp and its jacobian.
    TODO : set properly the list of arguments for this function, when
    things will be clearer ...
 */
typedef void (*ptrFunctionMCP)(int size , double* z, double * F);
typedef void (*ptrFunctionMCP2)(void* env, int n, double* z, double * F);
typedef void (*ptrFunctionMCP_nabla)(void* env, int n, double* z, NumericsMatrix * F);

/** \struct  MixedComplementarityProblem MixedComplementarityProblem.h
 * The structure that defines a Mixed Complementarity problem (MCP) : Find two vectors \f$(z,w \in {{\mathrm{I\!R}}}^{n+m})\f$ such that:\n
  \f{align*}{
  w &= \begin{pmatrix}w_e\\w_i\end{pmatrix} = F(z) \\
  w_e &=0 \\
  0 &\le w_i \perp z_i \ge 0
  \f}
  where "i" (resp. "e") stands for inequalities (resp. equalities). The vector \f$z\f$ is splitted like \f$w\f$:
  \f{equation*}{z =\begin{pmatrix}z_e\\z_i\end{pmatrix}\f}.
  \f$z_i,w_i\f$ are vectors of size <pre>sizeEqualities</pre>, \f$z_e,w_e\f$ vectors of size <pre>sizeInequalities</pre>
  and \f$F\f$ is a non linear function that must be user-defined.
 */
struct MixedComplementarityProblem
{
  int sizeEqualities; /**< size of equalities $z_e, w_e$ size */
  int sizeInequalities; /**< size of inequalities $z_i,w_i$ size */
  ptrFunctionMCP computeFmcp; /**< pointer to the function to compute F(z) */
  ptrFunctionMCP computeNablaFmcp; /** pointer to the function to compute the jacobian of F(z) */
  /** The value F(z) */
  double* Fmcp ;
  /** jacobian of F(z) */
  double* nablaFmcp ;

};

/** \struct MixedComplementarityProblem2 MixedComplementarityProblem.h
 * Structure that contains and defines a MixedComplementarityProblem
 */
struct MixedComplementarityProblem2
{
  int n1; /**< number of equalities constraints */
  int n2; /**< size of complementary variables */
  ptrFunctionMCP2 compute_Fmcp; /**< pointer to the function used to compute \f$F_{mcp}(z) = (G(z), H(z))\f$ */
  ptrFunctionMCP_nabla compute_nabla_Fmcp; /**< pointer to the function used to compute \f$\nabla_z F_{mcp}\f$ */
  NumericsMatrix* nabla_Fmcp; /**< storage for \f$\nabla_z F_{mcp}\f$*/
  void* env; /**< environment for the compute_Fmcp and compute_nabla_Fmcp function.
               When called from Python, it contains an object with compute_Fmcp and compute_nabla_Fmcp as methods.
               When called from C, it can reference a data struct containing variables needed for the computations.*/
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** \fn  void freeMixedComplementarityProblem(MixedComplementarityProblem* problem)
   *  \brief function to delete a MixedComplementarityProblem
   *  \param problem  pointer to a MixedComplementarityProblem to delete
   */
  void freeMixedComplementarityProblem(MixedComplementarityProblem* problem);

  /** free an MCP problem 
   * \param mcp structure to free
   */
  void freeMCP(MixedComplementarityProblem2* mcp);

  /** create an empty MCP problem
   * \return an MixedComplementarityProblem instance
   */
  MixedComplementarityProblem2* newMCP(void);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
