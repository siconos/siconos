// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef PLASTICITY2DPROBLEM_H
#define PLASTICITY2DPROBLEM_H

/*!\file Plasticity2DProblem.h
 * \brief Backward compatibility header - redirects to PlasticityProblem.h
 *
 * This header is kept for backward compatibility. New code should include
 * PlasticityProblem.h directly.
 *
 * The Plasticity2DProblem structure has been renamed to PlasticityProblem
 * and the eta/theta parameters have been moved to a model-specific structure
 * (Plasticity_DruckerPrager_model) accessible via problem->model.
 */

#include "PlasticityProblem.h"

#endif /* PLASTICITY2DPROBLEM_H */
