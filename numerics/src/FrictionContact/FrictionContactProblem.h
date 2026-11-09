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
#ifndef FRICTIONCONTACTPROBLEM_H
#define FRICTIONCONTACTPROBLEM_H

/*!\file FrictionContactProblem.h
  Definition of a structure to handle (reduced or dual) friction-contact
  (2D or 3D) problems.

  A friction-contact problem is defined by the data \f$(M, q, \mu)\f$.
  It consists in finding a local reaction \f$r\in{\mathrm{I\!R}}^m\f$ and a
  local relative velocity \f$u\in{\mathrm{I\!R}}^m\f$ such that
  \f[
    u = M r + q,
  \f]
  and, for each contact \f$\alpha=1,\dots,n_c\f$, the pair
  \f$(u_\alpha, r_\alpha)\f$ satisfies the Signorini--Coulomb friction law in
  dimension \f$d\f$ (with \f$d=2\f$ or \f$d=3\f$ and \f$m = d\,n_c\f$).

  This is the reduced (or dual) form obtained from a
  GlobalFrictionContactProblem by eliminating the global velocity, i.e.
  \f$M = H^T M_{\text{glob}}^{-1} H\f$ and
  \f$q = H^T M_{\text{glob}}^{-1} q_{\text{glob}} + b\f$.
*/
#include <stdio.h>  // for FILE

#include "NumericsFwd.h"     // for FrictionContactProblem, NumericsMatrix
#include "NumericsMatrix.h"  // for RawNumericsMatrix

/**
    The structure that defines a (reduced or dual) Friction-Contact (3D or 2D)
    problem.

    Given the data \f$(M, q, \mu)\f$, the unknowns are a local reaction
    \f$r\in{\mathrm{I\!R}}^m\f$ and a local relative velocity
    \f$u\in{\mathrm{I\!R}}^m\f$ satisfying
    \f[
      u = M r + q,
    \f]
    together with the Signorini--Coulomb friction law contact-wise in
    dimension \f$d\f$ (\f$d=2\f$ or \f$d=3\f$, \f$m = d\,n_c\f$).

    This reduced problem is obtained from a GlobalFrictionContactProblem by
    eliminating the global velocity:
    \f[
      M = H^T M_{\text{glob}}^{-1} H,\qquad
      q = H^T M_{\text{glob}}^{-1} q_{\text{glob}} + b.
    \f]
*/
struct FrictionContactProblem {
  /** dimension of the contact space (3D or 2D ) */
  int dimension;
  /** the number of contacts \f$ n_c \f$ */
  int numberOfContacts;
  /** \f$ {M} \in {{\mathrm{I\!R}}}^{m \times m} \f$,
     a matrix with \f$ m = d  n_c \f$ stored in NumericsMatrix structure */
  RawNumericsMatrix* M;
  /** \f$ {q} \in {{\mathrm{I\!R}}}^{m} \f$ */
  double* q;
  /** \f$ {\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of friction coefficients
      (\f$ n_c = \f$ numberOfContacts) */
  double* mu;
};

struct SplittedFrictionContactProblem {
  FrictionContactProblem* fc3d;
  NumericsMatrix* M_nn;
  NumericsMatrix* M_tn;
  NumericsMatrix* M_nt;
  NumericsMatrix* M_tt;
  double* q_n;
  double* q_t;
};

#if defined(__cplusplus)
extern "C" {
#endif

/* create an empty FrictionContactProblem
 * \return an empty fcp */
FrictionContactProblem* frictionContactProblem_new(void);

/** new FrictionContactProblem from minimal set of data
 *
 *  \param[in] dim the problem dimension
 *  \param[in] nc the number of contact
 *  \param[in] M the NumericsMatrix
 *  \param[in] q the q vector
 *  \param[in] mu the mu vector
 *  \return a pointer to a FrictionContactProblem structure
 */
FrictionContactProblem* frictionContactProblem_new_with_data(int dim, size_t nc,
                                                             NumericsMatrix* M, double* q,
                                                             double* mu);

/** free a FrictionContactProblem
 *
 *  \param problem the problem to free
 */
void frictionContactProblem_free(FrictionContactProblem* problem);

/** display a FrictionContactProblem
 *
 *  \param problem the problem to display
 */
void frictionContactProblem_display(FrictionContactProblem* problem);

/** print a FrictionContactProblem in a file (numerics .dat format)
 *
 *  \param problem the problem to print out
 *  \param file the dest file
 *  \return 0 if successfull
 */
int frictionContactProblem_printInFile(FrictionContactProblem* problem, FILE* file);

/** print a FrictionContactProblem in a file (numerics dat format)
 *
 *  \param problem the problem to print out
 *  \param filename the dest file
 *  \return 0 if successfull
 */
int frictionContactProblem_printInFilename(FrictionContactProblem* problem, char* filename);

/** read a FrictionContactProblem from a file descriptor
 *
 *  \param file descriptor
 *  \return problem the problem to read
 */
FrictionContactProblem* frictionContact_newFromFile(FILE* file);

/** read a FrictionContactProblem from a file (.dat or hdf5 if fclib is on) from
 *  its filename
 *
 *  \param filename the name of the input file
 *  \return problem the problem to read
 */
FrictionContactProblem* frictionContact_new_from_filename(const char* filename);

/** @brief create a SplittedFrictionContactProblem
 *
 *  @param problem the source problem
 *  @param splitted_problem the new (splitted) problem
 *  @return error code
 */
int createSplittedFrictionContactProblem(FrictionContactProblem* problem,
                                         SplittedFrictionContactProblem* splitted_problem);

/** @brief free a SplittedFrictionContactProblem
 *
 *  @param problem the problem to free
 */
void splittedFrictionContactProblem_free(SplittedFrictionContactProblem* problem);

void frictionContactProblem_compute_statistics(FrictionContactProblem* problem,
                                               double* reaction, double* velocity, double tol,
                                               int do_print);

/**
    Creates a new FrictionContact problem and initialize its content by copying
    an existing problem.

    \param problem the source problem to be copied
    \return a pointer to a new FrictionContactProblem
*/
FrictionContactProblem* frictionContact_copy(FrictionContactProblem* problem);

/**
    Rescales M matrix and q vector of a given FrictionContactProblem.

    \f[
    :math:`M = \alpha\gamma^2 M, q=\alpha\gamma q`
    \f]

    \param problem to be rescaled
    \param alpha rescaling factor
    \param gamma rescaling factor
*/
void frictionContact_rescaling(FrictionContactProblem* problem, double alpha, double gamma);

/** Status codes returned by FrictionContactProblem_compute_contact_status_3d_bnrt. */
enum SiconosFrictionContactStatus {
  SICONOS_FRICTION_CONTACT_STATUS_B =
      0, /**< Both velocity and reaction are in the friction cone. */
  SICONOS_FRICTION_CONTACT_STATUS_N, /**< Reaction is in the cone, velocity in the dual cone
                                        (no contact). */
  SICONOS_FRICTION_CONTACT_STATUS_R, /**< Sliding / rolling mode. */
  SICONOS_FRICTION_CONTACT_STATUS_T  /**< Take-off: sum is not in the friction cone. */
};

/** Classify a single 3D contact from its local velocity and reaction using the
 *  BNRT convention used by the IPM solvers.
 *
 *  B = both velocity and reaction are in the friction cone,
 *  N = reaction in cone, velocity in dual cone (non-contact),
 *  R = sliding/rolling regime,
 *  T = take-off (sum not in the cone).
 *
 *  \param[in] velocity local relative velocity of the contact (3 components)
 *  \param[in] reaction local reaction of the contact (3 components)
 *  \return a SiconosFrictionContactStatus value
 */
int FrictionContactProblem_compute_contact_status_3d_bnrt(const double velocity[3],
                                                          const double reaction[3]);

/** Fill an array with the BNRT status of each 3D contact.
 *
 *  \param[in] velocity local relative velocity array (3*n components)
 *  \param[in] reaction local reaction array (3*n components)
 *  \param[in] n number of contacts
 *  \param[out] statuses array of n ints, each a SiconosFrictionContactStatus value
 */
void FrictionContactProblem_get_contact_statuses_3d_bnrt(const double* velocity,
                                                         const double* reaction, size_t n,
                                                         int* statuses);

/** Count BNRT statuses for a set of 3D contacts.
 *
 *  \param[in] velocity local relative velocity array (3*n components)
 *  \param[in] reaction local reaction array (3*n components)
 *  \param[in] n number of contacts
 *  \param[out] nB,nN,nR,nT status counters (any may be NULL)
 */
void FrictionContactProblem_count_contact_statuses_3d_bnrt(const double* velocity,
                                                           const double* reaction, size_t n,
                                                           int* nB, int* nN, int* nR, int* nT);

/** Classify a 3D BNRT solution and print non-strict-complementarity cones.
 *
 *  \param[in] velocity local relative velocity array (3*n components)
 *  \param[in] reaction local reaction array (3*n components)
 *  \param[in] n number of contacts
 *  \param[out] nB,nN,nR BNRT counters (any may be NULL)
 *  \param[out] ns if all contacts are strict complementarity, the gap
 *                 (u+r)_0 - |(u+r)_bar| of the last contact
 *  \param[out] norm_sum if all contacts are strict complementarity, the
 *                 Euclidean norm of the last contact's u+r vector
 *  \param[in] verbose if nonzero, print detailed cone diagnostics for
 *                 non-strict-complementarity contacts
 *  \return the number of strict-complementarity contacts
 */
int FrictionContactProblem_compute_strict_complementarity_3d_bnrt(
    const double* velocity, const double* reaction, size_t n, int* nB, int* nN, int* nR,
    double* ns, double* norm_sum, int verbose);

/** Classify a single 2D contact from its local velocity, reaction and friction
 *  coefficient.
 *
 *  \return 0 = no-contact, 1 = stick, 2 = slip in the positive tangent
 *          direction, 3 = slip in the negative tangent direction.
 */
int FrictionContactProblem_compute_contact_status_2d(const double velocity[2],
                                                     const double reaction[2], double mu);

/** Classify a single 3D contact from its local velocity, reaction and friction
 *  coefficient.
 *
 *  \return 0 = no-contact, 1 = stick, 2 = slip.
 */
int FrictionContactProblem_compute_contact_status_3d(const double velocity[3],
                                                     const double reaction[3], double mu);

#if defined(__cplusplus)
}
#endif

#endif
