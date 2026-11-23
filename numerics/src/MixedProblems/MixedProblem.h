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

#ifndef MIXEDPROBLEM_H
#define MIXEDPROBLEM_H

/*!\file MixedProblem.h
  Definition of a structure to handle mixed problems combining
  friction-contact (2D) and plasticity models.
*/

#include <stdio.h>

#include "MixedProblems_options.h"
#include "NumericsFwd.h"
#include "NumericsMatrix.h"
#include "SiconosConfig.h"

/**
 * Enum to identify the type of model for each block in the mixed problem.
 */
enum MixedModelType {
  MIXED_MODEL_UNKNOWN = 0,
  MIXED_MODEL_FRICTION_2D,
  MIXED_MODEL_FRICTION_3D,
  MIXED_MODEL_ROLLING_FRICTION_2D,
  MIXED_MODEL_ROLLING_FRICTION_3D,
  MIXED_MODEL_VON_MISES,
  MIXED_MODEL_DRUCKER_PRAGER
};
typedef enum MixedModelType MixedModelType;

/**
 * Structure for Friction-Contact 2D model parameters.
 */
struct Mixed_Friction2D_model {
  double* mu;
};
typedef struct Mixed_Friction2D_model Mixed_Friction2D_model;

/**
 * Structure for Friction-Contact 3D model parameters.
 */
struct Mixed_Friction3D_model {
  double* mu;
};
typedef struct Mixed_Friction3D_model Mixed_Friction3D_model;

/**
 * Structure for Rolling Friction-Contact 2D model parameters.
 */
struct Mixed_RollingFriction2D_model {
  double* mu;  /* friction coefficient */
  double* mur; /* rolling friction coefficient */
};
typedef struct Mixed_RollingFriction2D_model Mixed_RollingFriction2D_model;

/**
 * Structure for Rolling Friction-Contact 3D model parameters.
 */
struct Mixed_RollingFriction3D_model {
  double* mu;  /* friction coefficient */
  double* mur; /* rolling friction coefficient */
};
typedef struct Mixed_RollingFriction3D_model Mixed_RollingFriction3D_model;

/**
 * Structure for Von Mises plasticity model parameters.
 */
struct Mixed_VonMises_model {
  double* sigma_y;
};
typedef struct Mixed_VonMises_model Mixed_VonMises_model;

/**
 * Structure for Drucker-Prager plasticity model parameters (see
 * Plasticity_DruckerPrager_model): K = { ||sigma_T|| <= eta sigma_N + c } for each cone, theta
 * the dilatancy coefficient. The arrays are not owned by the MixedProblem.
 */
struct Mixed_DruckerPrager_model {
  double* eta;
  double* theta;
  double* cohesion; /**< cohesion c per cone, NULL: no cohesion */
};
typedef struct Mixed_DruckerPrager_model Mixed_DruckerPrager_model;

/**
 * Union to hold model-specific parameters.
 */
union MixedModelUnion {
  Mixed_Friction2D_model* friction2d;
  Mixed_Friction3D_model* friction3d;
  Mixed_RollingFriction2D_model* rolling_friction2d;
  Mixed_RollingFriction3D_model* rolling_friction3d;
  Mixed_VonMises_model* von_mises;
  Mixed_DruckerPrager_model* drucker_prager;
  void* generic;
};
typedef union MixedModelUnion MixedModelUnion;

/**
 * Structure describing a single block/sub-problem in the mixed problem.
 * Each block has its own model type and parameters.
 */
struct MixedProblemBlock {
  /** Start index in the global vector (0-based) */
  int start_index;

  /** Dimension of this block (2 for friction 2D, 3 for plasticity 2D) */
  int dimension;

  /** Number of elements (contacts/cones) in this block */
  int num_elements;

  /** Model type for this block */
  MixedModelType model_type;

  /** Model-specific parameters */
  MixedModelUnion model;
};
typedef struct MixedProblemBlock MixedProblemBlock;

/**
 * The main structure that defines a Mixed problem.
 *
 * This structure combines multiple physical models (friction-contact and plasticity)
 * sharing a common matrix M and vector q.
 */
struct MixedProblem {
  /** Total dimension of the problem (sum of all block dimensions) */
  int total_dimension;

  /** Number of blocks/sub-problems */
  int num_blocks;

  /** Array of blocks describing each sub-problem */
  MixedProblemBlock* blocks;

  /**
   * Matrix M shared by all blocks.
   * Size: total_dimension x total_dimension
   */
  RawNumericsMatrix* M;

  /**
   * Vector q shared by all blocks.
   * Size: total_dimension
   */
  double* q;
};
typedef struct MixedProblem MixedProblem;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

/**
 * Create an empty MixedProblem.
 * \return a pointer to a new empty MixedProblem
 */
MixedProblem* mixedProblem_new(void);

/**
 * Create a MixedProblem with given data.
 *
 * \param[in] num_blocks number of blocks
 * \param[in] M the NumericsMatrix (shared)
 * \param[in] q the q vector (shared)
 * \return a pointer to a new MixedProblem
 */
MixedProblem* mixedProblem_new_with_data(int num_blocks, NumericsMatrix* M, double* q);

/**
 * Free a MixedProblem.
 *
 * \param[in,out] problem the problem to free
 */
void mixedProblem_free(MixedProblem* problem);

/**
 * Display a MixedProblem.
 *
 * \param[in] problem the problem to display
 */
void mixedProblem_display(MixedProblem* problem);

/**
 * Add a block to the mixed problem.
 *
 * \param[in,out] problem the mixed problem
 * \param[in] block_index index of the block to set (0 to num_blocks-1)
 * \param[in] start_index start index in global vector
 * \param[in] dimension dimension of this block
 * \param[in] num_elements number of elements
 * \param[in] model_type type of model
 * \return 0 if successful
 */
int mixedProblem_set_block(MixedProblem* problem, int block_index, int start_index,
                           int dimension, int num_elements, MixedModelType model_type);

/**
 * Set Friction-Contact 2D parameters for a block.
 *
 * \param[in,out] problem the mixed problem
 * \param[in] block_index index of the block
 * \param[in] mu array of friction coefficients (size: num_elements)
 * \return 0 if successful
 */
int mixedProblem_set_friction2d_params(MixedProblem* problem, int block_index, double* mu);

/**
 * Set Friction-Contact 3D parameters for a block.
 *
 * \param[in,out] problem the mixed problem
 * \param[in] block_index index of the block
 * \param[in] mu array of friction coefficients (size: num_elements)
 * \return 0 if successful
 */
int mixedProblem_set_friction3d_params(MixedProblem* problem, int block_index, double* mu);

/**
 * Set Rolling Friction-Contact 2D parameters for a block.
 *
 * \param[in,out] problem the mixed problem
 * \param[in] block_index index of the block
 * \param[in] mu array of friction coefficients (size: num_elements)
 * \param[in] mur array of rolling friction coefficients (size: num_elements)
 * \return 0 if successful
 */
int mixedProblem_set_rolling_friction2d_params(MixedProblem* problem, int block_index,
                                               double* mu, double* mur);

/**
 * Set Rolling Friction-Contact 3D parameters for a block.
 *
 * \param[in,out] problem the mixed problem
 * \param[in] block_index index of the block
 * \param[in] mu array of friction coefficients (size: num_elements)
 * \param[in] mur array of rolling friction coefficients (size: num_elements)
 * \return 0 if successful
 */
int mixedProblem_set_rolling_friction3d_params(MixedProblem* problem, int block_index,
                                               double* mu, double* mur);

/**
 * Set Von Mises parameters for a block.
 *
 * \param[in,out] problem the mixed problem
 * \param[in] block_index index of the block
 * \param[in] sigma_y array of yield stresses (size: num_elements)
 * \return 0 if successful
 */
int mixedProblem_set_vonmises_params(MixedProblem* problem, int block_index, double* sigma_y);

/**
 * Set Drucker-Prager parameters for a block.
 *
 * \param[in,out] problem the mixed problem
 * \param[in] block_index index of the block
 * \param[in] eta array of cone coefficients (size: num_elements)
 * \param[in] theta array of dilatancy coefficients (size: num_elements)
 * \return 0 if successful
 */
int mixedProblem_set_druckerprager_params(MixedProblem* problem, int block_index, double* eta,
                                          double* theta);

/**
 * Set the cohesion of a Drucker-Prager block (after mixedProblem_set_druckerprager_params).
 * The cone of element e is then { ||sigma_T|| <= eta_e sigma_N + c_e } (eta_e > 0 if
 * c_e != 0).
 *
 * \param[in,out] problem the mixed problem
 * \param[in] block_index index of the block
 * \param[in] cohesion array of cohesions (size: num_elements, not owned), NULL: no cohesion
 * \return 0 if successful
 */
int mixedProblem_set_druckerprager_cohesion(MixedProblem* problem, int block_index,
                                            double* cohesion);

/**
 * Get the string representation of a model type.
 *
 * \param[in] model_type the model type
 * \return string description
 */
const char* mixedModelType_to_string(MixedModelType model_type);

/**
 * General driver for Mixed Problems.
 *
 * \param[in] problem the MixedProblem to solve
 * \param[in,out] z the solution vector (size: problem->total_dimension)
 * \param[in,out] w the auxiliary vector (size: problem->total_dimension)
 * \param[in,out] options solver options structure
 * \return 0 if successful, error code otherwise
 */
int mixedProblem_driver(MixedProblem* problem, double* z, double* w, SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif /* MIXEDPROBLEM_H */
