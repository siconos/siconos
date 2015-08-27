/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef SOCLCP_H
#define SOCLCP_H

/*! \page soclcpProblem
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
 *    C^\alpha_{\mu} = \{ r \mid \|[r_1,\ldot,r_{n^\alpha}]\| \leq \mu^\alpha * r_0   \} \subset {\mathrm{I\!R}}}^{n^\alpha}
 * \f}
 *
 *  \section SOCLCPSolversList Available solvers for SOCCLP
 * Use the generic function SecondOrderConeLinearComplementarityProblem_driver() to call one the the specific solvers listed below:
 *
 * <ul>
 *
 * <li> SecondOrderConeLinearComplementarityProblem_psor(): PSOR (Gauss-Seidel with overrelaxation) solver.
 *       SolverId : SICONOS_SOCLCP_NSGS = , </li>
 *
 *
 * </ul>
 * (see the functions/solvers list in SecondOrderConeLinearComplementarityProblem_Solvers.h)
 *
 *
 */


/*!\file SecondOrderConeLinearComplementarityProblem.h
  \brief Definition of a structure to handle with the SecondOrderConeLinearComplementarityProblem
*/

#include "NumericsMatrix.h"

/** \struct  SecondOrderConeLinearComplementarityProblem SecondOrderConeLinearComplementarityProblem.h
 *  The structure that defines a Second Order Cone Linear Complementarity Problem
 *  \f$\mathrm{SOCLCP}(M,q,\mu)\f$  such that
 * \f{eqnarray*}{
 * \begin{cases}
 *   u = M r + q \\
 *   C_{\mu} = \{ r \} =  \prod_{\alpha =1}^{n_c} C^\alpha_{\mu} \\
 *   C^\alpha_{\mu} = \{ r \mid \|[r_1,\ldot,r_{n^\alpha}]\| \leq \mu^\alpha * r_0   \} \subset {\mathrm{I\!R}}}^{n^\alpha}
 *   C^\star_{\mu} \ni  {u} \perp r \in C_{\mu}
 * \end{cases}
 * \f}
 *   \param nc the number of contacts \f$ n_c \f$
 *   \param M \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
 *    a matrix with \f$ n = d  n_c\f$ stored in NumericsMatrix structure
 *   \param q  \f${q} \in {{\mathrm{I\!R}}}^{n} \f$,
 *   coneDimensions \f${\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of dimension of the cones
 *     (\f$ n_c =\f$ numberOfContacts)
 *   \param mu \f${\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of  coefficients
 *      (\f$ n_c =\f$ numberOfContacts)
 */
typedef struct
{
  /** the number of cones \f$ n_c \f$ in the Cartesian product */
  int nc;
  /** M \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
     a matrix with \f$ n = d  n_c\f$ stored in NumericsMatrix structure */
  NumericsMatrix* M;
  /** \f${q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double* q;
  /** coneDimensions \f${\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of dimension of the cones
      (\f$ n_c =\f$ numberOfContacts) */
  int* coneDimensions;
  /** mu \f${\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of coefficients
      (\f$ n_c =\f$ nc) */
  double* mu;
} SecondOrderConeLinearComplementarityProblem;


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
   * \param[in] nc the number of contact
   * \param[in] M the NumericsMatrix
   * \param[in] q the q vector
   * \param[in] coneDimensions the q vector
   * \param[in] mu the mu vector
   * \return a pointer to a SecondOrderConeLinearComplementarityProblem structure
   */
  SecondOrderConeLinearComplementarityProblem* secondOrderConeLinearComplementarityProblem_new( int nc,
                                                                          NumericsMatrix* M, double* q,
                                                                          int *coneDimensions, double* mu);



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
