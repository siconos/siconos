// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef NM_MA57_h
#define NM_MA57_h

#if defined(__cplusplus)
extern "C" {
#endif

#ifdef WITH_MA57
#include <lbl.h>

/** Free the config data for MA57
 * \param param p a pointer on the linear solver parameters
 */
void NM_MA57_free(void* p);

#endif /* WITH_MA57 */

#if defined(__cplusplus)
}
#endif

#endif
