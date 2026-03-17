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
