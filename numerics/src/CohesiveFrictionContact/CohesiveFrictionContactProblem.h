/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef COHESIVEFRICTIONCONTACTPROBLEM_H
#define COHESIVEFRICTIONCONTACTPROBLEM_H

/*!\file CohesiveFrictionContactProblem.h
  Definition of a structure to handle cohesive friction-contact (2D or 3D) problems.

  The cohesive friction-contact problem extends the standard friction-contact
  problem by adding cohesive forces that act before contact occurs. This models
  the behavior of interfaces that can sustain traction before failure.

  The problem is formulated as:
  - velocity = q + M * reaction + V * cohesion
  - friction law with cohesion-modified gap
  - cohesion forces computed from internal state variables

  \see FrictionContactProblem for the base problem structure
  \see RollingFrictionContactProblem for a similar extension pattern
*/

#include <stdio.h>  // for FILE

#include "NumericsFwd.h"     // for CohesiveFrictionContactProblem, NumericsMatrix
#include "NumericsMatrix.h"  // for RawNumericsMatrix

/**
    The structure that defines a Cohesive Friction-Contact (3D or 2D) problem.

    This extends the standard FrictionContactProblem by adding cohesion parameters
    that model the attractive forces between surfaces before contact occurs.
*/
struct CohesiveFrictionContactProblem {
  /** dimension of the contact space (3D = 3 or 2D = 2) */
  int dimension;
  /** the number of contacts \f$ n_c \f$ */
  int numberOfContacts;
  /** Number of cohesive points \f$ n_{coh} \f$ (may differ from numberOfContacts) */
  int numberOfCohesivePoints;
  
  /** \f$ {M} \in {{
      \mathrm{I\!R}}}^{m \times m} \f$,
     a matrix with \f$ m = d  (n_c + n_{coh}} \f$ stored in NumericsMatrix structure */
  RawNumericsMatrix *M;
  
  /** \f$ {q} \in {{
      \mathrm{I\!R}}}^{m} \f$ - the local velocity and displacement vector */
  double *q;
  
  /** \f$ {\mu} \in {{
      \mathrm{I\!R}}}^{n_c} \f$, vector of friction coefficients
      (\f$ n_c =\f$ numberOfContacts) */
  double *mu;

  /** \f$ {c_n} \in {{
      \mathrm{I\!R}}}^{n_{coh}} \f$, vector of cohesion intensity in normal direction */
  double *c_n;
  /** \f$ {c_t} \in {{
      \mathrm{I\!R}}}^{n_{coh}} \f$, vector of cohesion intensity in tangent direction */
  double *c_t;
  
  /** \f$ {q_v} \in {{
      \mathrm{I\!R}}}^{n} \f with \f$ m = d  n \f$ $,  vector associated with cohesive points
  */
  double *q_v;

  /** \f$ {q_u} \in {{
      \mathrm{I\!R}}}^{n} \f with \f$ m = d  n_{coh} \f$ $,  vector associated with cohesive points
  */
  double *q_u;
  
  /** Matrix W for mapping cohesive forces to contact space (required).
      \f$ {V} \in {{
      \mathrm{I\!R}}}^{m \times m} \f$*/
  RawNumericsMatrix *W;
   
  /** Matrix W for mapping cohesive forces to contact space (required).
      \f$ {V} \in {{
      \mathrm{I\!R}}}^{m \times n} \f$*/
  
  RawNumericsMatrix *V;
  /** Matrix X for additional cohesive contributions (required).
      Used for extended cohesive zone models.
      \f$ {X} \in {{
      \mathrm{I\!R}}}^{n \times n} \f$*/
  RawNumericsMatrix *X;

  /** Matrix U for coupling terms (required).
      Used for coupling between different contact modes.
      \f$ {X} \in {{
      \mathrm{I\!R}}}^{n \times m} \f$*/
  RawNumericsMatrix *U;
  
  /** Internal state variables for cohesive zone models (optional).
      This stores damage parameters, displacement history, etc. */
  double *internal_state;
  /** Size of internal state vector per contact */
  int internal_state_size;
};

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * Create an empty CohesiveFrictionContactProblem
 * \return an empty cohesive friction-contact problem
 */
CohesiveFrictionContactProblem *cohesiveFrictionContactProblem_new(void);

/** new CohesiveFrictionContactProblem from minimal set of data
 *
 *  \param[in] dim the problem dimension (2 or 3)
 *  \param[in] nc the number of contact
 *  \param[in] n_coh the number of cohesive points
 *  \param[in] M the NumericsMatrix
 *  \param[in] q the q vector
 *  \param[in] mu the mu vector (friction coefficients)
 *  \param[in] V the V matrix (required)
 *  \param[in] X the X matrix (required)
 *  \param[in] U the U matrix (required)
 *  \return a pointer to a CohesiveFrictionContactProblem structure
 */
CohesiveFrictionContactProblem *cohesiveFrictionContactProblem_new_with_data(
    int dim, int nc, int n_coh, NumericsMatrix *M, double *q, double *mu, NumericsMatrix *V,
    NumericsMatrix *X, NumericsMatrix *U);

/** free a CohesiveFrictionContactProblem
 *
 *  \param problem the problem to free
 */
void cohesiveFrictionContactProblem_free(CohesiveFrictionContactProblem *problem);

/** display a CohesiveFrictionContactProblem
 *
 *  \param problem the problem to display
 */
void cohesiveFrictionContact_display(CohesiveFrictionContactProblem *problem);

/** print a CohesiveFrictionContactProblem in a file (numerics .dat format)
 *
 *  \param problem the problem to print out
 *  \param file the destination file
 *  \return 0 if successful
 */
int cohesiveFrictionContact_printInFile(CohesiveFrictionContactProblem *problem, FILE *file);

/** print a CohesiveFrictionContactProblem in a file (numerics .dat format) from
 *  its filename
 *
 *  \param problem the problem to print out
 *  \param filename the destination file
 *  \return 0 if successful
 */
int cohesiveFrictionContact_printInFilename(CohesiveFrictionContactProblem *problem,
                                            char *filename);

/** read a CohesiveFrictionContactProblem from a file (numerics .dat format)
 *
 *  \param file the source file
 *  \return a pointer to a CohesiveFrictionContactProblem structure
 */
CohesiveFrictionContactProblem *cohesiveFrictionContact_newFromFile(FILE *file);

/** read a CohesiveFrictionContactProblem from a file (numerics .dat format) from
 *  its filename
 *
 *  \param filename the source file
 *  \return a pointer to a CohesiveFrictionContactProblem structure
 */
CohesiveFrictionContactProblem *cohesiveFrictionContact_newFromFilename(char *filename);

/**
 * Update the cohesion intensity vectors in the problem
 *
 * This allows updating the cohesive intensities based on the current state
 * of the simulation (e.g., after each time step).
 *
 * \param problem the cohesive friction contact problem
 * \param c_n the new normal cohesion intensity vector
 * \param c_t the new tangent cohesion intensity vector
 */
void cohesiveFrictionContactProblem_set_cohesion(CohesiveFrictionContactProblem *problem,
                                                 double *c_n, double *c_t);

/**
 * Compute the effective q vector including cohesive contribution
 *
 * q_eff = q + V * c (where c combines c_n and c_t)
 *
 * \param problem the cohesive friction contact problem
 * \param q_eff output vector for effective q (must be pre-allocated)
 */
void cohesiveFrictionContactProblem_compute_effective_q(
    CohesiveFrictionContactProblem *problem, double *q_eff);

/**
 * Build the global matrix M and vector q from block components
 *
 * Constructs the block matrix M = [[W, V], [U, X]] and vector q = [q_v, q_u]
 * from the individual block components using NM_insert.
 *
 * \param[in,out] problem the cohesive friction contact problem
 * \return 0 if successful, error code otherwise
 *
 * \note This function allocates M and q if they are NULL, or reuses them if already allocated.
 *       The matrices W, V, U, X and vectors q_v, q_u must be set in the problem before calling.
 */
int cohesiveFrictionContactProblem_build_M_q_from_blocks(
    CohesiveFrictionContactProblem *problem);

#if defined(__cplusplus)
}
#endif

#endif  // COHESIVEFRICTIONCONTACTPROBLEM_H
