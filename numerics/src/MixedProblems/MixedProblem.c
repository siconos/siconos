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

#include "MixedProblem.h"

#include <stdlib.h>
#include <string.h>

#include "MixedProblem_NSGS.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"

MixedProblem* mixedProblem_new(void) {
  MixedProblem* problem = (MixedProblem*)malloc(sizeof(MixedProblem));
  if (problem == NULL) {
    int err = numerics_error("mixedProblem_new", "Failed to allocate MixedProblem");
    (void)err;  // Error already printed
    return NULL;
  }

  problem->total_dimension = 0;
  problem->num_blocks = 0;
  problem->blocks = NULL;
  problem->M = NULL;
  problem->q = NULL;

  return problem;
}

MixedProblem* mixedProblem_new_with_data(int num_blocks, NumericsMatrix* M, double* q) {
  MixedProblem* problem = mixedProblem_new();
  if (problem == NULL) return NULL;

  problem->num_blocks = num_blocks;
  problem->M = M;
  problem->q = q;

  if (num_blocks > 0) {
    problem->blocks = (MixedProblemBlock*)calloc(num_blocks, sizeof(MixedProblemBlock));
    if (problem->blocks == NULL) {
      int err =
          numerics_error("mixedProblem_new_with_data", "Failed to allocate blocks array");
      (void)err;  // Error already printed, just need to clean up
      free(problem);
      return NULL;
    }

    /* Initialize all blocks to unknown type */
    for (int i = 0; i < num_blocks; i++) {
      problem->blocks[i].model_type = MIXED_MODEL_UNKNOWN;
      problem->blocks[i].start_index = 0;
      problem->blocks[i].dimension = 0;
      problem->blocks[i].num_elements = 0;
      problem->blocks[i].model.generic = NULL;
    }
  }

  return problem;
}

void mixedProblem_free(MixedProblem* problem) {
  if (problem == NULL) return;

  if (problem->blocks != NULL) {
    /* Free model-specific data for each block */
    for (int i = 0; i < problem->num_blocks; i++) {
      MixedProblemBlock* block = &problem->blocks[i];

      switch (block->model_type) {
        case MIXED_MODEL_FRICTION_2D: {
          if (block->model.friction2d != NULL) {
            free(block->model.friction2d);
          }
          break;
        }
        case MIXED_MODEL_FRICTION_3D: {
          if (block->model.friction3d != NULL) {
            free(block->model.friction3d);
          }
          break;
        }
        case MIXED_MODEL_ROLLING_FRICTION_2D: {
          if (block->model.rolling_friction2d != NULL) {
            free(block->model.rolling_friction2d);
          }
          break;
        }
        case MIXED_MODEL_ROLLING_FRICTION_3D: {
          if (block->model.rolling_friction3d != NULL) {
            free(block->model.rolling_friction3d);
          }
          break;
        }
        case MIXED_MODEL_VON_MISES: {
          if (block->model.von_mises != NULL) {
            free(block->model.von_mises);
          }
          break;
        }
        case MIXED_MODEL_DRUCKER_PRAGER: {
          if (block->model.drucker_prager != NULL) {
            free(block->model.drucker_prager);
          }
          break;
        }
        default:
          break;
      }
    }

    free(problem->blocks);
  }

  /* Note: M and q are not freed here as they are typically shared/externally managed */

  free(problem);
}

void mixedProblem_display(MixedProblem* problem) {
  if (problem == NULL) {
    printf("MixedProblem: NULL\n");
    return;
  }

  printf("MixedProblem:\n");
  printf("  Total dimension: %d\n", problem->total_dimension);
  printf("  Number of blocks: %d\n", problem->num_blocks);
  printf("  Matrix M: %p\n", (void*)problem->M);
  printf("  Vector q: %p\n", (void*)problem->q);
  printf("  Blocks:\n");

  for (int i = 0; i < problem->num_blocks; i++) {
    MixedProblemBlock* block = &problem->blocks[i];
    printf("    Block %d:\n", i);
    printf("      Start index: %d\n", block->start_index);
    printf("      Dimension: %d\n", block->dimension);
    printf("      Num elements: %d\n", block->num_elements);
    printf("      Model type: %s\n", mixedModelType_to_string(block->model_type));

    switch (block->model_type) {
      case MIXED_MODEL_FRICTION_2D: {
        if (block->model.friction2d != NULL && block->model.friction2d->mu != NULL) {
          printf("      Friction coefficients (mu):\n");
          for (int j = 0; j < block->num_elements; j++) {
            printf("        [%d]: %e\n", j, block->model.friction2d->mu[j]);
          }
        }
        break;
      }
      case MIXED_MODEL_FRICTION_3D: {
        if (block->model.friction3d != NULL && block->model.friction3d->mu != NULL) {
          printf("      Friction coefficients (mu):\n");
          for (int j = 0; j < block->num_elements; j++) {
            printf("        [%d]: %e\n", j, block->model.friction3d->mu[j]);
          }
        }
        break;
      }
      case MIXED_MODEL_ROLLING_FRICTION_2D: {
        if (block->model.rolling_friction2d != NULL) {
          if (block->model.rolling_friction2d->mu != NULL) {
            printf("      Friction coefficients (mu):\n");
            for (int j = 0; j < block->num_elements; j++) {
              printf("        [%d]: %e\n", j, block->model.rolling_friction2d->mu[j]);
            }
          }
          if (block->model.rolling_friction2d->mur != NULL) {
            printf("      Rolling friction coefficients (mur):\n");
            for (int j = 0; j < block->num_elements; j++) {
              printf("        [%d]: %e\n", j, block->model.rolling_friction2d->mur[j]);
            }
          }
        }
        break;
      }
      case MIXED_MODEL_ROLLING_FRICTION_3D: {
        if (block->model.rolling_friction3d != NULL) {
          if (block->model.rolling_friction3d->mu != NULL) {
            printf("      Friction coefficients (mu):\n");
            for (int j = 0; j < block->num_elements; j++) {
              printf("        [%d]: %e\n", j, block->model.rolling_friction3d->mu[j]);
            }
          }
          if (block->model.rolling_friction3d->mur != NULL) {
            printf("      Rolling friction coefficients (mur):\n");
            for (int j = 0; j < block->num_elements; j++) {
              printf("        [%d]: %e\n", j, block->model.rolling_friction3d->mur[j]);
            }
          }
        }
        break;
      }
      case MIXED_MODEL_VON_MISES: {
        if (block->model.von_mises != NULL && block->model.von_mises->sigma_y != NULL) {
          printf("      Yield stresses (sigma_y):\n");
          for (int j = 0; j < block->num_elements; j++) {
            printf("        [%d]: %e\n", j, block->model.von_mises->sigma_y[j]);
          }
        }
        break;
      }
      case MIXED_MODEL_DRUCKER_PRAGER: {
        if (block->model.drucker_prager != NULL) {
          if (block->model.drucker_prager->eta != NULL) {
            printf("      Cone coefficients (eta):\n");
            for (int j = 0; j < block->num_elements; j++) {
              printf("        [%d]: %e\n", j, block->model.drucker_prager->eta[j]);
            }
          }
          if (block->model.drucker_prager->theta != NULL) {
            printf("      Dilatancy coefficients (theta):\n");
            for (int j = 0; j < block->num_elements; j++) {
              printf("        [%d]: %e\n", j, block->model.drucker_prager->theta[j]);
            }
          }
        }
        break;
      }
      default:
        break;
    }
  }
}

int mixedProblem_set_block(MixedProblem* problem, int block_index, int start_index,
                           int dimension, int num_elements, MixedModelType model_type) {
  if (problem == NULL) {
    return numerics_error("mixedProblem_set_block", "NULL problem pointer");
  }

  if (block_index < 0 || block_index >= problem->num_blocks) {
    return numerics_error("mixedProblem_set_block", "Invalid block index");
  }

  MixedProblemBlock* block = &problem->blocks[block_index];

  block->start_index = start_index;
  block->dimension = dimension;
  block->num_elements = num_elements;
  block->model_type = model_type;
  block->model.generic = NULL;

  /* Update total dimension */
  if (start_index + dimension * num_elements > problem->total_dimension) {
    problem->total_dimension = start_index + dimension * num_elements;
  }

  return 0;
}

int mixedProblem_set_friction2d_params(MixedProblem* problem, int block_index, double* mu) {
  if (problem == NULL) {
    return numerics_error("mixedProblem_set_friction2d_params", "NULL problem pointer");
  }

  if (block_index < 0 || block_index >= problem->num_blocks) {
    return numerics_error("mixedProblem_set_friction2d_params", "Invalid block index");
  }

  MixedProblemBlock* block = &problem->blocks[block_index];

  if (block->model_type != MIXED_MODEL_FRICTION_2D) {
    return numerics_error("mixedProblem_set_friction2d_params",
                          "Block is not a Friction-Contact 2D model");
  }

  Mixed_Friction2D_model* model =
      (Mixed_Friction2D_model*)malloc(sizeof(Mixed_Friction2D_model));
  if (model == NULL) {
    return numerics_error("mixedProblem_set_friction2d_params",
                          "Failed to allocate Friction2D model");
  }

  model->mu = mu;
  block->model.friction2d = model;

  return 0;
}

int mixedProblem_set_friction3d_params(MixedProblem* problem, int block_index, double* mu) {
  if (problem == NULL) {
    return numerics_error("mixedProblem_set_friction3d_params", "NULL problem pointer");
  }

  if (block_index < 0 || block_index >= problem->num_blocks) {
    return numerics_error("mixedProblem_set_friction3d_params", "Invalid block index");
  }

  MixedProblemBlock* block = &problem->blocks[block_index];

  if (block->model_type != MIXED_MODEL_FRICTION_3D) {
    return numerics_error("mixedProblem_set_friction3d_params",
                          "Block is not a Friction-Contact 3D model");
  }

  Mixed_Friction3D_model* model =
      (Mixed_Friction3D_model*)malloc(sizeof(Mixed_Friction3D_model));
  if (model == NULL) {
    return numerics_error("mixedProblem_set_friction3d_params",
                          "Failed to allocate Friction3D model");
  }

  model->mu = mu;
  block->model.friction3d = model;

  return 0;
}

int mixedProblem_set_rolling_friction2d_params(MixedProblem* problem, int block_index,
                                               double* mu, double* mur) {
  if (problem == NULL) {
    return numerics_error("mixedProblem_set_rolling_friction2d_params",
                          "NULL problem pointer");
  }

  if (block_index < 0 || block_index >= problem->num_blocks) {
    return numerics_error("mixedProblem_set_rolling_friction2d_params", "Invalid block index");
  }

  MixedProblemBlock* block = &problem->blocks[block_index];

  if (block->model_type != MIXED_MODEL_ROLLING_FRICTION_2D) {
    return numerics_error("mixedProblem_set_rolling_friction2d_params",
                          "Block is not a Rolling Friction-Contact 2D model");
  }

  Mixed_RollingFriction2D_model* model =
      (Mixed_RollingFriction2D_model*)malloc(sizeof(Mixed_RollingFriction2D_model));
  if (model == NULL) {
    return numerics_error("mixedProblem_set_rolling_friction2d_params",
                          "Failed to allocate RollingFriction2D model");
  }

  model->mu = mu;
  model->mur = mur;
  block->model.rolling_friction2d = model;

  return 0;
}

int mixedProblem_set_rolling_friction3d_params(MixedProblem* problem, int block_index,
                                               double* mu, double* mur) {
  if (problem == NULL) {
    return numerics_error("mixedProblem_set_rolling_friction3d_params",
                          "NULL problem pointer");
  }

  if (block_index < 0 || block_index >= problem->num_blocks) {
    return numerics_error("mixedProblem_set_rolling_friction3d_params", "Invalid block index");
  }

  MixedProblemBlock* block = &problem->blocks[block_index];

  if (block->model_type != MIXED_MODEL_ROLLING_FRICTION_3D) {
    return numerics_error("mixedProblem_set_rolling_friction3d_params",
                          "Block is not a Rolling Friction-Contact 3D model");
  }

  Mixed_RollingFriction3D_model* model =
      (Mixed_RollingFriction3D_model*)malloc(sizeof(Mixed_RollingFriction3D_model));
  if (model == NULL) {
    return numerics_error("mixedProblem_set_rolling_friction3d_params",
                          "Failed to allocate RollingFriction3D model");
  }

  model->mu = mu;
  model->mur = mur;
  block->model.rolling_friction3d = model;

  return 0;
}

int mixedProblem_set_vonmises_params(MixedProblem* problem, int block_index, double* sigma_y) {
  if (problem == NULL) {
    return numerics_error("mixedProblem_set_vonmises_params", "NULL problem pointer");
  }

  if (block_index < 0 || block_index >= problem->num_blocks) {
    return numerics_error("mixedProblem_set_vonmises_params", "Invalid block index");
  }

  MixedProblemBlock* block = &problem->blocks[block_index];

  if (block->model_type != MIXED_MODEL_VON_MISES) {
    return numerics_error("mixedProblem_set_vonmises_params",
                          "Block is not a Von Mises model");
  }

  Mixed_VonMises_model* model = (Mixed_VonMises_model*)malloc(sizeof(Mixed_VonMises_model));
  if (model == NULL) {
    return numerics_error("mixedProblem_set_vonmises_params",
                          "Failed to allocate Von Mises model");
  }

  model->sigma_y = sigma_y;
  block->model.von_mises = model;

  return 0;
}

int mixedProblem_set_druckerprager_params(MixedProblem* problem, int block_index, double* eta,
                                          double* theta) {
  if (problem == NULL) {
    return numerics_error("mixedProblem_set_druckerprager_params", "NULL problem pointer");
  }

  if (block_index < 0 || block_index >= problem->num_blocks) {
    return numerics_error("mixedProblem_set_druckerprager_params", "Invalid block index");
  }

  MixedProblemBlock* block = &problem->blocks[block_index];

  if (block->model_type != MIXED_MODEL_DRUCKER_PRAGER) {
    return numerics_error("mixedProblem_set_druckerprager_params",
                          "Block is not a Drucker-Prager model");
  }

  Mixed_DruckerPrager_model* model =
      (Mixed_DruckerPrager_model*)malloc(sizeof(Mixed_DruckerPrager_model));
  if (model == NULL) {
    return numerics_error("mixedProblem_set_druckerprager_params",
                          "Failed to allocate Drucker-Prager model");
  }

  model->eta = eta;
  model->theta = theta;
  model->cohesion = NULL;
  block->model.drucker_prager = model;

  return 0;
}

int mixedProblem_set_druckerprager_cohesion(MixedProblem* problem, int block_index,
                                            double* cohesion) {
  if (problem == NULL) {
    return numerics_error("mixedProblem_set_druckerprager_cohesion", "NULL problem pointer");
  }
  if (block_index < 0 || block_index >= problem->num_blocks) {
    return numerics_error("mixedProblem_set_druckerprager_cohesion", "Invalid block index");
  }
  MixedProblemBlock* block = &problem->blocks[block_index];
  if (block->model_type != MIXED_MODEL_DRUCKER_PRAGER || !block->model.drucker_prager) {
    return numerics_error("mixedProblem_set_druckerprager_cohesion",
                          "Block is not a Drucker-Prager model with parameters set");
  }
  if (cohesion) {
    for (int e = 0; e < block->num_elements; e++)
      if (cohesion[e] != 0.0 && !(block->model.drucker_prager->eta[e] > 0.0))
        return numerics_error("mixedProblem_set_druckerprager_cohesion",
                              "element %d: cohesion with eta <= 0", e);
  }
  block->model.drucker_prager->cohesion = cohesion;
  return 0;
}

const char* mixedModelType_to_string(MixedModelType model_type) {
  switch (model_type) {
    case MIXED_MODEL_UNKNOWN:
      return "UNKNOWN";
    case MIXED_MODEL_FRICTION_2D:
      return "FRICTION_2D";
    case MIXED_MODEL_FRICTION_3D:
      return "FRICTION_3D";
    case MIXED_MODEL_ROLLING_FRICTION_2D:
      return "ROLLING_FRICTION_2D";
    case MIXED_MODEL_ROLLING_FRICTION_3D:
      return "ROLLING_FRICTION_3D";
    case MIXED_MODEL_VON_MISES:
      return "VON_MISES";
    case MIXED_MODEL_DRUCKER_PRAGER:
      return "DRUCKER_PRAGER";
    default:
      return "INVALID";
  }
}

/* ===========================================================================
 * Driver for Mixed Problems
 * =========================================================================== */

int mixedProblem_driver(MixedProblem* problem, double* z, double* w, SolverOptions* options) {
  if (problem == NULL || z == NULL || w == NULL || options == NULL) {
    return numerics_error("mixedProblem_driver", "NULL pointer argument");
  }

  int info = NUMERICS_OK;

  /* Initialize based on solver ID */
  switch (options->solverId) {
    case SICONOS_MIXED_NSGS:
      mixedProblem_NSGS_init(problem, options);
      if (options->solverData == NULL) {
        return NUMERICS_ERR_MEMORY;
      }
      mixedProblem_NSGS(problem, z, w, &info, options);
      mixedProblem_NSGS_free(problem, options);
      break;

    default: {
      int err = numerics_error("mixedProblem_driver", "Unknown solver ID");
      (void)err;
    }
      return NUMERICS_ERR_INVALID_SOLVER;
  }

  return info;
}
