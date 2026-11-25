// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file vertex_extraction.h
 * \brief interface to various LP solvers
 */

#ifndef vertex_extraction_h
#define vertex_extraction_h

#include "SiconosLapack.h"  // IWYU pragma: keep - For lapack_int
#include "SiconosSets.h"    // for polyhedron

#if defined(__cplusplus)
extern "C" {
#endif

void siconos_find_vertex(const polyhedron* P, unsigned size, lapack_int* basis);

#if defined(__cplusplus)
}
#endif

#endif
