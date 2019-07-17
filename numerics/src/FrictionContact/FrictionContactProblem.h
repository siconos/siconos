/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
*/

#include "NumericsFwd.h"
#include "SiconosConfig.h"
#include "NumericsMatrix.h"

/** \struct FrictionContactProblem FrictionContactProblem.h
 *  The structure that defines a (reduced or dual) Friction-Contact (3D or 2D) problem.
*/
struct FrictionContactProblem
{
  /** dimension of the contact space (3D or 2D ) */
  int dimension;
  /** the number of contacts \f$ n_c \f$ */
  int numberOfContacts;
  /** \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
     a matrix with \f$ n = d  n_c\f$ stored in NumericsMatrix structure */
  RawNumericsMatrix* M;
  /** \f${q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double* q;
  /** \f${\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of friction coefficients
      (\f$ n_c =\f$ numberOfContacts) */
  double* mu;
};

struct SplittedFrictionContactProblem
{
  FrictionContactProblem * fc3d;
  NumericsMatrix * M_nn;
  NumericsMatrix * M_tn;
  NumericsMatrix * M_nt;
  NumericsMatrix * M_tt;
  double * q_n;
  double * q_t;
};




#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif


  /* create an empty FrictionContactProblem
   * \return an empty fcp */
  FrictionContactProblem* frictionContactProblem_new(void);

  /** new FrictionContactProblem from minimal set of data
   * \param[in] dim the problem dimension
   * \param[in] nc the number of contact
   * \param[in] M the NumericsMatrix
   * \param[in] q the q vector
   * \param[in] mu the mu vector
   * \return a pointer to a FrictionContactProblem structure
   */
  FrictionContactProblem* frictionContactProblem_new_with_data(int dim, int nc,
      NumericsMatrix* M, double* q, double* mu);

  /** free a FrictionContactProblem
   * \param problem the problem to free
   */
  void frictionContactProblem_free(FrictionContactProblem* problem);


  /** display a FrictionContactProblem
   * \param problem the problem to display
   */
  void frictionContact_display(FrictionContactProblem*  problem);

  /** print a FrictionContactProblem in a file (numerics .dat format)
   * \param problem the problem to print out
   * \param file the dest file
   * \return 0 if successfull
   */
  int frictionContact_printInFile(FrictionContactProblem*  problem, FILE* file);

  /** print a FrictionContactProblem in a file (numerics .dat format) from its filename
   * \param problem the problem to print out
   * \param filename the dest file
   * \return 0 if successfull
   */
  int frictionContact_printInFilename(FrictionContactProblem*  problem, char * filename);

  /** read a FrictionContactProblem in a file (numerics .dat format)
   * \param problem the problem to read
   * \param file the target file
   * \return 0 if successfull
   */
  int frictionContact_newFromFile(FrictionContactProblem*  problem, FILE* file);

  /** read a FrictionContactProblem in a file (numerics .dat format) from its filename
   * \param problem the problem to read
   * \param filename the name of the target file
   * \return 0 if successfull
   */
  int frictionContact_newFromFilename(FrictionContactProblem*  problem, char * filename);






  void createSplittedFrictionContactProblem(FrictionContactProblem* problem, SplittedFrictionContactProblem * splitted_problem);

  void frictionContactProblem_compute_statistics(FrictionContactProblem* problem,
                                                 double * reaction,
                                                 double * velocity,
                                                 double tol,
                                                 int do_print)  ;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
