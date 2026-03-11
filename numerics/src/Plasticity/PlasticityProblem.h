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
#ifndef PLASTICITYPROBLEM_H
#define PLASTICITYPROBLEM_H

/*!\file PlasticityProblem.h
  Definition of a structure to handle Plasticity problems.
*/
#include <stdio.h>  // for FILE

#include "NumericsFwd.h"     // for NumericsMatrix
#include "NumericsMatrix.h"  // for RawNumericsMatrix
#include "SiconosConfig.h"   // for BUILD_AS_CPP // IWYU pragma: keep

/**
    Drucker-Prager plasticity model parameters.
    This structure contains the model-specific parameters for Drucker-Prager
    plasticity (which includes Mohr-Coulomb as a special case).
*/
struct Plasticity_DruckerPrager_model {
  /** \f$ {\eta} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of cone coefficients
      (\f$ n_c = \f$ numberOfCones) */
  double *eta;
  /** \f$ {\theta} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of dilatancy coefficients
      (\f$ n_c = \f$ numberOfCones) */
  double *theta;
};

typedef struct Plasticity_DruckerPrager_model Plasticity_DruckerPrager_model;

/**
    The structure that defines a Plasticity problem.
    This is a generic structure that can be used with different plasticity models.
*/
struct PlasticityProblem {
  /** dimension of the stress space */
  int dimension;
  /** the number of  \f$ n_c \f$ */
  int numberOfCones;
  /** \f$ {M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
     a matrix with \f$ n = d  n_c \f$ stored in NumericsMatrix structure */
  RawNumericsMatrix *M;
  /** \f$ {q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double *q;
  /** Pointer to model-specific parameters (e.g., Drucker-Prager model).
      NULL if no model is set. */
  Plasticity_DruckerPrager_model *model;
};

typedef struct PlasticityProblem PlasticityProblem;

/* Backward compatibility: Plasticity2DProblem is now PlasticityProblem */
typedef PlasticityProblem Plasticity2DProblem;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/**
    Create a new Drucker-Prager model with given parameters.
    
    \param[in] eta vector of cone coefficients
    \param[in] theta vector of dilatancy coefficients
    \return a pointer to a new Plasticity_DruckerPrager_model structure
*/
Plasticity_DruckerPrager_model *plasticity_DruckerPrager_model_new(double *eta, double *theta);

/**
    Free a Drucker-Prager model.
    
    \param[in] model the model to free
*/
void plasticity_DruckerPrager_model_free(Plasticity_DruckerPrager_model *model);

/* create an empty PlasticityProblem
 * \return an empty problem */
PlasticityProblem *plasticityProblem_new(void);

/** new PlasticityProblem from minimal set of data
 *
 *  \param[in] dim the problem dimension
 *  \param[in] nc the number of contact
 *  \param[in] M the NumericsMatrix
 *  \param[in] q the q vector
 *  \param[in] model pointer to model parameters (can be NULL)
 *  \return a pointer to a PlasticityProblem structure
 */
PlasticityProblem *plasticityProblem_new_with_data(int dim, int nc, NumericsMatrix *M,
                                                    double *q, Plasticity_DruckerPrager_model *model);

/** free a PlasticityProblem
 *
 *  \param problem the problem to free
 */
void plasticityProblem_free(PlasticityProblem *problem);

/** display a PlasticityProblem
 *
 *  \param problem the problem to display
 */
void plasticity_display(PlasticityProblem *problem);

/** print a PlasticityProblem in a file (numerics .dat format)
 *
 *  \param problem the problem to print out
 *  \param file the dest file
 *  \return 0 if successfull
 */
int plasticity_printInFile(PlasticityProblem *problem, FILE *file);

/** print a PlasticityProblem in a file (numerics dat format)
 *
 *  \param problem the problem to print out
 *  \param filename the dest file
 *  \return 0 if successfull
 */
int plasticity_printInFilename(PlasticityProblem *problem, char *filename);

/** read a PlasticityProblem from a file descriptor
 *
 *  \param file descriptor
 *  \return problem the problem to read
 */
PlasticityProblem *plasticity_newFromFile(FILE *file);

/** read a PlasticityProblem from a file (.dat or hdf5 if fclib is on) from
 *  its filename
 *
 *  \param filename the name of the input file
 *  \return problem the problem to read
 */
PlasticityProblem *plasticity_new_from_filename(const char *filename);

/**
    Creates a new Plasticity problem and initialize its content by copying
    an existing problem.

    \param problem the source problem to be copied
    \return a pointer to a new PlasticityProblem
*/
PlasticityProblem *plasticity_copy(PlasticityProblem *problem);

/**
    Rescales M matrix and q vector of a given PlasticityProblem.

    \f[
    :math:`M = \alpha\gamma^2 M, q=\alpha\gamma q`
    \f]

    \param problem to be rescaled
    \param alpha rescaling factor
    \param gamma rescaling factor
*/
void plasticity_rescaling(PlasticityProblem *problem, double alpha, double gamma);

/* Backward compatibility function names (deprecated) */
Plasticity2DProblem *plasticity2DProblem_new(void);
Plasticity2DProblem *plasticity2DProblem_new_with_data(int dim, int nc, NumericsMatrix *M,
                                                         double *q, double *eta, double *theta);
void plasticity2DProblem_free(Plasticity2DProblem *problem);
void plasticity2D_display(Plasticity2DProblem *problem);
int plasticity2D_printInFile(Plasticity2DProblem *problem, FILE *file);
int plasticity2D_printInFilename(Plasticity2DProblem *problem, char *filename);
Plasticity2DProblem *plasticity2D_newFromFile(FILE *file);
Plasticity2DProblem *plasticity2D_new_from_filename(const char *filename);
Plasticity2DProblem *plasticity2D_copy(Plasticity2DProblem *problem);
void plasticity2D_rescaling(Plasticity2DProblem *problem, double alpha, double gamma);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
