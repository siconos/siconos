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
#ifndef MOHRCOULOMB2DPROBLEM_H
#define MOHRCOULOMB2DPROBLEM_H

/*!\file MohrCoulomb2DProblem.h
  Definition of a structure to handle Mohr Coulomb 2D plasticity problems.
*/
#include <stdio.h>  // for FILE

#include "NumericsFwd.h"     // for NumericsMatrix
#include "NumericsMatrix.h"  // for RawNumericsMatrix
#include "SiconosConfig.h"   // for BUILD_AS_CPP // IWYU pragma: keep

/**
    The structure that defines a Mohr Coulomb 2D xproblem.
*/
struct MohrCoulomb2DProblem {
  /** dimension of the stress space */
  int dimension;
  /** the number of  \f$ n_c \f$ */
  int numberOfCones;
  /** \f$ {M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
     a matrix with \f$ n = d  n_c \f$ stored in NumericsMatrix structure */
  RawNumericsMatrix *M;
  /** \f$ {q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double *q;
  /** \f$ {\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of cone coefficients
      (\f$ n_c = \f$ numberOfCones) */
  double *eta;
  /** \f$ {\theta} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of dilatency coefficients
      (\f$ n_c = \f$ numberOfCones) */
  double *theta;
};

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/* create an empty MohrCoulomb2DProblem
 * \return an empty fcp */
MohrCoulomb2DProblem *mohrCoulomb2DProblem_new(void);

/** new MohrCoulomb2DProblem from minimal set of data
 *
 *  \param[in] dim the problem dimension
 *  \param[in] nc the number of contact
 *  \param[in] M the NumericsMatrix
 *  \param[in] q the q vector
 *  \param[in] theta the theta vector
 *  \return a pointer to a MohrCoulomb2DProblem structure
 */
MohrCoulomb2DProblem *mohrCoulomb2DProblem_new_with_data(int dim, int nc, NumericsMatrix *M,
                                                         double *q, double *eta, double *theta);

/** free a MohrCoulomb2DProblem
 *
 *  \param problem the problem to free
 */
void mohrCoulomb2DProblem_free(MohrCoulomb2DProblem *problem);

/** display a MohrCoulomb2DProblem
 *
 *  \param problem the problem to display
 */
void mohrCoulomb2D_display(MohrCoulomb2DProblem *problem);

/** print a MohrCoulomb2DProblem in a file (numerics .dat format)
 *
 *  \param problem the problem to print out
 *  \param file the dest file
 *  \return 0 if successfull
 */
int mohrCoulomb2D_printInFile(MohrCoulomb2DProblem *problem, FILE *file);

/** print a MohrCoulomb2DProblem in a file (numerics dat format)
 *
 *  \param problem the problem to print out
 *  \param filename the dest file
 *  \return 0 if successfull
 */
int mohrCoulomb2D_printInFilename(MohrCoulomb2DProblem *problem, char *filename);

/** read a MohrCoulomb2DProblem from a file descriptor
 *
 *  \param file descriptor
 *  \return problem the problem to read
 */
MohrCoulomb2DProblem *mohrCoulomb2D_newFromFile(FILE *file);

/** read a MohrCoulomb2DProblem from a file (.dat or hdf5 if fclib is on) from
 *  its filename
 *
 *  \param filename the name of the input file
 *  \return problem the problem to read
 */
MohrCoulomb2DProblem *mohrCoulomb2D_new_from_filename(const char *filename);

/**
    Creates a new MohrCoulomb2D problem and initialize its content by copying
    an existing problem.

    \param problem the source problem to be copied
    \return a pointer to a new MohrCoulomb2DProblem
*/
MohrCoulomb2DProblem *mohrCoulomb2D_copy(MohrCoulomb2DProblem *problem);

/**
    Rescales M matrix and q vector of a given MohrCoulomb2DProblem.

    \f[
    :math:`M = \alpha\gamma^2 M, q=\alpha\gamma q`
    \f]

    \param problem to be rescaled
    \param alpha rescaling factor
    \param gamma rescaling factor
*/
void mohrCoulomb2D_rescaling(MohrCoulomb2DProblem *problem, double alpha, double gamma);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
