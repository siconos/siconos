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

/*! \file GenericMechanicalProblem.h
 * \brief GenericMechanicalProblem and GMP_LocalProblem structures.
 */

#ifndef NUMERICSGENERICMECHANICALPROBLEM_H
#define NUMERICSGENERICMECHANICALPROBLEM_H

#include <stdio.h>

#include "NumericsFwd.h"

/* void * solverFC3D; */
/* void * solverEquality; */
/* void * solverLCP; */
/* void * solverMLCP; */

/** \struct GMP_LocalProblem GenericMechanicalProblem.h
 *  \brief One local sub-problem inside a GenericMechanicalProblem.
 *
 *  A GenericMechanicalProblem is a list of local problems coupled through a
 *  global block matrix M. Each GMP_LocalProblem stores the type-specific
 *  formulation (LCP, equality, friction-contact, ...) and its own local
 *  right-hand side q_local.
 *
 *  \param type       problem type (see SICONOS_NUMERICS_PROBLEM_TYPE)
 *  \param problem    type-specific problem struct (e.g. LinearComplementarityProblem)
 *  \param q_local    local right-hand side, owned by this struct
 *  \param size       dimension of the local problem
 *  \param error      non-zero if the local solver reported an error
 *  \param next       next local problem in the list
 *  \param previous   previous local problem in the list
 */
struct GMP_LocalProblem {
  int type;
  void* problem;
  double* q_local; /* local right-hand side, owned by this struct */
  size_t size;     /* size of the local problem */
  int error;       /* non-zero if the local solver reported an error */
  struct GMP_LocalProblem* next;
  struct GMP_LocalProblem* previous;
};

/** Backward-compatible alias for GMP_LocalProblem. */
typedef struct GMP_LocalProblem listNumericsProblem;

/** \enum SICONOS_NUMERICS_PROBLEM_TYPE ids for the possible/allowed numerics problem
 * formulations
 */
typedef enum {
  SICONOS_NUMERICS_PROBLEM_LCP = 0,
  SICONOS_NUMERICS_PROBLEM_MLCP = 1,
  SICONOS_NUMERICS_PROBLEM_EQUALITY = 2,
  SICONOS_NUMERICS_PROBLEM_FC2D = 3,
  SICONOS_NUMERICS_PROBLEM_FC3D = 4,
  SICONOS_NUMERICS_PROBLEM_NCP = 5,
  SICONOS_NUMERICS_PROBLEM_MCP = 6,
  SICONOS_NUMERICS_PROBLEM_VI = 7,
  SICONOS_NUMERICS_PROBLEM_AVI = 8,
  SICONOS_NUMERICS_PROBLEM_RELAY = 9,
} SICONOS_NUMERICS_PROBLEM_TYPE;

enum NUMERICS_GMP_FREE { GMP_FREE_MATRIX = 4, GMP_FREE_GMP = 8 };

/** \struct GenericMechanicalProblem GenericMechanicalProblem.h
 * \brief A mixed non-smooth mechanical problem.
 *
 * \param globalSize    total size of the global problem (sum of local sizes)
 * \param maxLocalSize  maximal size of a local problem
 * \param M             global NumericsMatrix (set by the user)
 * \param q             global right-hand side vector (set by the user)
 * \param firstLocal    first local problem in the list (private, managed by gmp_add)
 * \param lastLocal     last local problem in the list (private, managed by gmp_add)
 *
 *  ONLY M and q must be allocated/freed by the user; the other fields are
 *  private. Do not fill this structure by hand: use genericMechanicalProblem_new(),
 *  gmp_add() and genericMechanicalProblem_free().
 */
struct GenericMechanicalProblem {
  /* Total size of the global problem (sum of all local sizes). */
  /* PRIVATE: managed by gmp_add. */
  size_t globalSize;
  /* Maximal size of a local problem. */
  /* PRIVATE: managed by gmp_add. */
  size_t maxLocalSize;
  /* Global matrix, must be set by the user. */
  NumericsMatrix* M;
  /* Global right-hand side vector, must be set by the user. */
  double* q;
  /* First local problem in the list. PRIVATE: managed by gmp_add. */
  GMP_LocalProblem* firstLocal;
  /* Last local problem in the list. PRIVATE: managed by gmp_add. */
  GMP_LocalProblem* lastLocal;
};

#if defined(__cplusplus)
extern "C" {
#endif

/** Build an empty GenericMechanicalProblem.
 * \return a pointer on the built GenericMechanicalProblem.
 */
GenericMechanicalProblem* genericMechanicalProblem_new(void);

/** Free the list of contained sub-problems and, depending on \a level, the
 *  global matrix and vector. Also frees \a problem.
 *
 *  \param[in,out] problem   the GenericMechanicalProblem to free
 *  \param[in]     level  GMP_FREE_GMP to free the problem only, or
 *                        GMP_FREE_MATRIX to also free M and q
 */
void genericMechanicalProblem_free(GenericMechanicalProblem* problem, unsigned int level);

/** Print a GenericMechanicalProblem to a file.
 *  \param[in] problem the printed problem
 *  \param[in,out] file the output file
 */
void genericMechanicalProblem_printInFile(GenericMechanicalProblem* problem, FILE* file);

/** Read a GenericMechanicalProblem from a file descriptor.
 * \param[in] file the file descriptor
 * \return the read problem, or NULL on error
 */
GenericMechanicalProblem* genericMechanical_newFromFile(FILE* file);

/** Read a GenericMechanicalProblem from a file name.
 * \param[in] filename the name of the input file
 * \return the read problem, or NULL on error
 */
GenericMechanicalProblem* genericMechanical_new_from_filename(const char* filename);

/** Display a GenericMechanicalProblem on standard output.
 *  \param[in] problem the displayed problem
 */
void genericMechanicalProblem_display(GenericMechanicalProblem* problem);

/** Insert a local problem into a GenericMechanicalProblem.
 *
 *  The memory of the elementary block matrix is not managed here: the user
 *  must ensure it. In Siconos, the Kernel allocates the global NumericsMatrix
 *  and the diagonal block is plugged later in gmp_gauss_seidel
 *  (localProblem->M->matrix0 = diagBlock).
 *
 *  \param[in,out] problem         the GenericMechanicalProblem
 *  \param[in]     problemType  type of the added sub-problem
 *  \param[in]     size         size of the local formulation
 *  \return a pointer to the type-specific problem struct, or NULL on error
 */
void* gmp_add(GenericMechanicalProblem* problem, SICONOS_NUMERICS_PROBLEM_TYPE problemType,
              size_t size);

/** Return the non-smooth problem formulation name from its id.
 *  \param[in] id problem id (must be a valid SICONOS_NUMERICS_PROBLEM_TYPE)
 *  \return the problem name, or NULL if unknown
 */
const char* ns_problem_id_to_name(SICONOS_NUMERICS_PROBLEM_TYPE id);

#if defined(__cplusplus)
}
#endif

#endif
