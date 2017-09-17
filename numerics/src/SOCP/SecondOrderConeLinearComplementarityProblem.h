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
#ifndef SOCLCP_H
#define SOCLCP_H

/*! \page soclcpProblem Second Order Cone Linear Complementarity Problem
 *
 * \section soclcpIntro Problem statement
 *  Given
 * <ul>
 *   <li> a symmetric positive semi--definite  matrix \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$ </li>
 *   <li> a vector \f$ {q} \in {{\mathrm{I\!R}}}^n\f$</li>
 *   <li> a vector of coefficients\f$\mu \in{{\mathrm{I\!R}}}^{n_c}\f$</li>
 *</ul>
 * the second order cone linear complementarity problem (SOCLCP) is to find two vectors \f$u\in{{\mathrm{I\!R}}}^n\f$,
 * and \f$r\in {{\mathrm{I\!R}}}^n\f$,
 * denoted by \f$\mathrm{SOCCLP}(M,q,\mu)\f$  such that
 * \f{eqnarray*}{
 * \begin{cases}
 *   u = M r + q \\                               \
 *    C^\star_{\mu} \ni {u} \perp r \in C_{\mu}
 * \end{cases}
 * \f}
 * and the set \f$C^{\alpha,\star}_{\mu^\alpha}\f$ is its dual.
 *
 * The set C is the second order cone given by
 * \f{eqnarray}{
 *    C_{\mu} = \{ r \} =  \prod_{\alpha =1}^{n_c} C^\alpha_{\mu}
 * \f}
 * with
 * \f{eqnarray}{
 *    C^\alpha_{\mu} = \{ r \mid \|[r_1, \ldots, r_{n^\alpha}]\| \leq \mu^\alpha * r_0   \} \subset {\mathrm{I\!R}}^{n^\alpha}
 * \f}
 *
 * The problem is stored and given to the solver in numerics thanks to
 *  the C structure  SecondOrderConeLinearComplementarityProblem.

 *  \section SOCLCPSolversList Available solvers for SOCCLP

see SOCLCP_cst.h for solver ids.
 
 * Use the generic function soclcp_driver() to call one the the specific solvers listed below:
 *
 * <ul>
 *
 * <li> soclcp_nsgs() : PSOR (Gauss-Seidel with overrelaxation) solver.
 *       SolverId : SICONOS_SOCLCP_NSGS , </li>
 * <li> soclcp_VI_FixedPointProjection() : VI formulation and fixed point projection.
 *       SolverId : SICONOS_SOCLCP_VI_FPP , </li>
 * <li> soclcp_VI_ExtraGradient() : VI formulation and extra-gradient solver.
 *       SolverId : SICONOS_SOCLCP_VI_EG , </li>
  * </ul>
 * See the related functions/solvers list in SOCLCP_Solvers.h.
 *
 */


/*!\file SecondOrderConeLinearComplementarityProblem.h
  \brief Definition of a structure to handle with the SecondOrderConeLinearComplementarityProblem
*/

#include "NumericsMatrix.h"
#include "NumericsFwd.h"

/** \struct  SecondOrderConeLinearComplementarityProblem SecondOrderConeLinearComplementarityProblem.h
 *  The structure that defines a Second Order Cone Linear Complementarity Problem, , see details in \ref soclcpProblem.
 */
struct SecondOrderConeLinearComplementarityProblem
{
  /** the problem dimension. must be equal to\f$ \sum_{i}^{n_c} d_i \f$   */
  int n;

  /** the number of cones \f$ n_c \f$ in the Cartesian product */
  int nc;
  /** \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
     a matrix with \f$ n = d  n_c\f$ stored in NumericsMatrix structure */
  NumericsMatrix* M;
  /** \f${q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double* q;
  /** \f${coneIndex} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of indices of the cones
      (\f$ n_c =\f$ nc) */
  unsigned int* coneIndex;
  /** \f${\tau} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of coefficients
      (\f$ n_c =\f$ nc) */
  double* tau;
};


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"

{
#endif
/** display a SecondOrderConeLinearComplementarityProblem
 * \param problem the problem to display
 */
void secondOrderConeLinearComplementarityProblem_display(SecondOrderConeLinearComplementarityProblem*  problem);

/** print a SecondOrderConeLinearComplementarityProblem in a file (numerics .dat format)
 * \param problem the problem to print out
 * \param file the dest file
 * \return 0 if successfull
 */
int secondOrderConeLinearComplementarityProblem_printInFile(SecondOrderConeLinearComplementarityProblem*  problem, FILE* file);

/** print a SecondOrderConeLinearComplementarityProblem in a file (numerics .dat format) from its filename
 * \param problem the problem to print out
 * \param filename the dest file
 * \return 0 if successfull
 */
int secondOrderConeLinearComplementarityProblem_printInFilename(SecondOrderConeLinearComplementarityProblem*  problem, char * filename);

/** read a SecondOrderConeLinearComplementarityProblem in a file (numerics .dat format)
 * \param problem the problem to read
 * \param file the target file
 * \return 0 if successfull
 */
int secondOrderConeLinearComplementarityProblem_newFromFile(SecondOrderConeLinearComplementarityProblem*  problem, FILE* file);

/** read a SecondOrderConeLinearComplementarityProblem in a file (numerics .dat format) from its filename
 * \param problem the problem to read
 * \param filename the name of the target file
 * \return 0 if successfull
 */
int secondOrderConeLinearComplementarityProblem_newFromFilename(SecondOrderConeLinearComplementarityProblem*  problem, char * filename);

/** free a SecondOrderConeLinearComplementarityProblem
 * \param problem the problem to free
 */
void freeSecondOrderConeLinearComplementarityProblem(SecondOrderConeLinearComplementarityProblem* problem);


/** new SecondOrderConeLinearComplementarityProblem from minimal set of data
 * \param[in] n the size of the problem
 * \param[in] nc the number of contact
 * \param[in] M the NumericsMatrix
 * \param[in] q the q vector
 * \param[in] coneIndex
 * \param[in] mu the mu vector
 * \return a pointer to a SecondOrderConeLinearComplementarityProblem structure
 */
SecondOrderConeLinearComplementarityProblem* secondOrderConeLinearComplementarityProblem_new
(int n, int nc, NumericsMatrix* M, double* q,
 unsigned int *coneIndex, double* mu);



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
