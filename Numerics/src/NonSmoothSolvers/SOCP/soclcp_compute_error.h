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
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef SOCLCP_compute_error_H
#define SOCLCP_compute_error_H

/*!\file soclcp_compute_error.h
  \brief functions related to error computation for SOCLCP

  \author Vincent Acary, 28/08/2015

*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

/** Error computation for SOCLCP problem
    \param problem the structure which defines the SOCLCP
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param[in,out] error value
    \return 0 if ok
 */
int soclcp_compute_error(SecondOrderConeLinearComplementarityProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);

/** Error computation for one SOCLCP problem
    \param z vector
    \param w vector
    \param mu coeficient of friction
    \param[in,out] error value
 */
void soclcp_unitary_compute_and_add_error(double z[3] , double w[3], double mu, double * error);

/** Error computation for SOCLCP problem
    \param problem the structure which defines the SOCLCP
    \param z vector
    \param w vector
    \param options
    \param tolerance value for error computation
    \param[in,out] error value
    \return 0 if ok
 */
int soclcp_compute_error_velocity(SecondOrderConeLinearComplementarityProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);



/** Error computation for SOCLCP problem with Tresca Friction
    \param problem the structure which defines the SOCLCP
    \param z vector
    \param w vector
    \param tolerance value for error computation
    \param options
    \param[in,out] error value
    \return 0 if ok
 */
int soclcp_Tresca_compute_error(SecondOrderConeLinearComplementarityProblem* problem, double *z , double *w, double tolerance, SolverOptions * options, double * error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
