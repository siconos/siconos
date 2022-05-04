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
#ifndef GLOBALROLLINGFRICTIONCONTACTPROBLEM_H
#define GLOBALROLLINGFRICTIONCONTACTPROBLEM_H

/*!\file GlobalRollingFrictionContactProblem.h
  Definition of a structure to handle friction-contact (2D or 3D)
  problems.
*/

#include "NumericsFwd.h" // for GlobalRollingFrictionContactProblem, NumericsMatrix
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#include <stdio.h> // for FILE

/** \struct GlobalRollingFrictionContactProblem
 * GlobalRollingFrictionContactProblem.h The structure that defines a (reduced
 * or dual) Friction-Contact (3D or 2D) problem.
 */
struct GlobalRollingFrictionContactProblem {
  /** dimension of the contact space (3D or 2D ) */
  int dimension;
  /** the number of contacts \f$ n_c \f$ */
  int numberOfContacts;
  /** \f$ {M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
     a matrix with \f$ n = d  n_c \f$ stored in NumericsMatrix structure */
  NumericsMatrix *M;
  /**  \f$ {H} \in {{\mathrm{I\!R}}}^{n \times m} \f$,
       a matrix with \f$ m = d  n_c \f$ stored in NumericsMatrix structure */
  NumericsMatrix *H;
  /** \f$ {q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double *q;
  /** \f$ {b} \in {{\mathrm{I\!R}}}^{m} \f$ */
  double *b;
  /** \f$ {\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of friction coefficients
      (\f$ n_c = \f$ numberOfContacts) */
  double *mu;
  /** \f$ {\mu_r} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of friction
     coefficients
      (\f$ n_c = \f$ numberOfContacts) */
  double *mu_r;
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/* create an empty GlobalRollingFrictionContactProblem
 * \return an empty fcp */
GlobalRollingFrictionContactProblem *
globalRollingFrictionContactProblem_new(void);

/** new GlobalRollingFrictionContactProblem from minimal set of data
 *
 *  \param[in] dim the problem dimension
 *  \param[in] nc the number of contact
 *  \param[in] M the NumericsMatrix
 *  \param[in] q the q vector
 *  \param[in] mu the mu vector
 *  \return a pointer to a GlobalRollingFrictionContactProblem structure
 */
GlobalRollingFrictionContactProblem *
globalRollingFrictionContactProblem_new_with_data(int dim, int nc,
                                                  NumericsMatrix *M, double *q,
                                                  double *mu, double *mu_r);

/** free a GlobalRollingFrictionContactProblem
 *
 *  \param problem the problem to free
 */
void globalRollingFrictionContactProblem_free(
    GlobalRollingFrictionContactProblem *problem);

/** display a GlobalRollingFrictionContactProblem
 *
 *  \param problem the problem to display
 */
void globalRollingFrictionContact_display(
    GlobalRollingFrictionContactProblem *problem);

/** print a GlobalRollingFrictionContactProblem in a file (numerics .dat format)
 *
 *  \param problem the problem to print out
 *  \param file the dest file
 *  \return 0 if successfull
 */
int globalRollingFrictionContact_printInFile(
    GlobalRollingFrictionContactProblem *problem, FILE *file);

/** print a GlobalRollingFrictionContactProblem in a file (numerics .dat format)
 *  from its filename
 *
 *  \param problem the problem to print out
 *  \param filename the dest file
 *  \return 0 if successfull
 */
int globalRollingFrictionContact_printInFilename(
    GlobalRollingFrictionContactProblem *problem, char *filename);

/** read a GlobalRollingFrictionContactProblem from a file descriptor
 *
 *  \param file descriptor
 *  \return problem the problem to read
 */
GlobalRollingFrictionContactProblem *
globalRollingFrictionContact_newFromFile(FILE *file);

/** read a GlobalRollingFrictionContactProblem from a file (.dat or hdf5 if
 *  fclib is on) from its filename
 *
 *  \param filename the name of the input file
 *  \return problem the problem to read
 */
GlobalRollingFrictionContactProblem *
globalRollingFrictionContact_new_from_filename(const char *filename);

/** Compute the global velocity given the reaction
 *
 *  \param[in] problem to be considered
 *  \param[in] reaction the reaction, if there is no contacts
 *  reaction can be NULL
 *  \param[out] globalVelocity the global velocity computed by inverting the
 *  system.
 */
int globalRollingFrictionContact_computeGlobalVelocity(
    GlobalRollingFrictionContactProblem *problem, double *reaction,
    double *globalVelocity);

/* Reformulation into local problem */
RollingFrictionContactProblem *
globalRollingFrictionContact_reformulation_RollingFrictionContact(
    GlobalRollingFrictionContactProblem *problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
