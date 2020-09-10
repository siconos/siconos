/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#ifndef GLOBALFRICTIONCONTACTPROBLEM_H
#define GLOBALFRICTIONCONTACTPROBLEM_H

/*!\file GlobalFrictionContactProblem.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
*/


#include <stdio.h>        // for NULL, FILE
#include "NumericsFwd.h"  // for GlobalFrictionContactProblem, NumericsMatrix
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep
#include "NumericsMatrix.h" // for BalancingMatrices

/** \struct GlobalFrictionContactProblem GlobalFrictionContactProblem.h
 *
 * The structure that defines a Friction-Contact (3D or 2D ) problem
 *
 \rst
 
 Details in :ref:`global_fc_problem`.
 \endrst

*/

struct GlobalFrictionContactProblem
{
  /** dimension \f$d=2\f$ or \f$d=3\f$ of the contact space (3D or 2D ) */
  int dimension;
  /** the number of contacts \f$ n_c \f$ */
  int numberOfContacts;
  /** \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
      a matrix with \f$ n\f$ stored in NumericsMatrix structure */
  NumericsMatrix* M;
  /**  \f${H} \in {{\mathrm{I\!R}}}^{n \times m} \f$,
       a matrix with \f$ m = d  n_c\f$ stored in NumericsMatrix structure */
  NumericsMatrix* H;
  /** \f${q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double* q;
  /** \f${b} \in {{\mathrm{I\!R}}}^{m} \f$ */
  double* b;
  /** mu \f${\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of friction coefficients
      (\f$ n_c =\f$ numberOfContacts) */
  double* mu;
  /** opaque environment, solver specific */
  void* env;
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /* creates an empty GlobalFrictionContactProblem
     \return a pointer to a GlobalFrictionContactProblem
  */
  GlobalFrictionContactProblem* globalFrictionContactProblem_new(void);

  /** displays the problem onto screen
   * \param[in] problem to be displayed
   */
  void globalFrictionContact_display(GlobalFrictionContactProblem*  problem);


  /** Saves problem struct into a file.
      \param[in] problem structure
      \param[in] file, file descriptor
      \return file status (1 if everything has worked properly)
  */
  int globalFrictionContact_printInFile(GlobalFrictionContactProblem*  problem, FILE* file);

  /** Saves problem struct into a file.
      \param[in] problem structure
      \param[in] filename, name of the input file
      \return file status (1 if everything has worked properly)
  */
  int globalFrictionContact_printInFileName(GlobalFrictionContactProblem*  problem,
                                            const char * filename);

  /** read a GlobalFrictionContactProblem from a file descriptor
   * \param file descriptor
   * \return problem the problem to read
   */
  GlobalFrictionContactProblem*  globalFrictionContact_newFromFile(FILE* file);

  /** read a GlobalFrictionContactProblem from a file (.dat or hdf5 if fclib is on) from its filename
   * \param filename the name of the input file
   * \return problem the problem to read
   */
  GlobalFrictionContactProblem* globalFrictionContact_new_from_filename(const char * filename);

  /** Release memory for the problem structure
      \param[inout] problem, global-Friction problem structure to be freed.
  */
  void globalFrictionContact_free(GlobalFrictionContactProblem* problem);

  GlobalFrictionContactProblem* globalFrictionContact_copy(GlobalFrictionContactProblem* problem);

  void globalFrictionContact_rescaling(GlobalFrictionContactProblem* problem, double alpha,  double beta, double gamma);
  void globalFrictionContact_balancing_M(
    GlobalFrictionContactProblem* problem,
    BalancingMatrices * B_for_M);
  void globalFrictionContact_balancing_M_H(
    GlobalFrictionContactProblem* problem,
    BalancingMatrices * B_for_M,
    BalancingMatrices * B_for_H);

  
  /** Compute the global velocity given the reaction
   * \param[in] problem to be considered
   * \param[in] reaction the reaction, if there is no contacts reaction can be NULL
   * \param[out] globalVelocity the global velocity computed by inverting the system.
   */
  int globalFrictionContact_computeGlobalVelocity(
    GlobalFrictionContactProblem* problem,
    double * reaction,
    double * globalVelocity);



#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif
#endif
