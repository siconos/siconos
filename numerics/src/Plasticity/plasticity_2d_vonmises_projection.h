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
#ifndef PLASTICITY_2D_VONMISES_PROJECTION_H
#define PLASTICITY_2D_VONMISES_PROJECTION_H

/*!\file plasticity_2d_vonmises_projection.h
 * \brief Projection functions for Von Mises plasticity model
 */

#include "PlasticityProblem.h"

#if defined(__cplusplus)
extern "C" {
#endif

/**
    Project a trial stress onto the Von Mises yield surface.
    
    The Von Mises yield criterion in 2D (plane strain) is:
    f(σ) = √(σ_x² + σ_y² - σ_x*σ_y + 3*τ_xy²) - σ_y ≤ 0
    
    This uses the radial return algorithm (closest point projection).
    
    \param[in,out] stress trial stress vector [σ_x, σ_y, τ_xy]
    \param[in] sigma_y yield stress
    \return 0 if projection successful, 1 if stress was inside yield surface (no projection needed)
*/
int plasticity_2d_projectionOnVonMises(double stress[3], double sigma_y);

/**
    Compute the Von Mises equivalent stress (q).
    
    q = √(σ_x² + σ_y² - σ_x*σ_y + 3*τ_xy²)
    
    \param[in] stress stress vector [σ_x, σ_y, τ_xy]
    \return Von Mises equivalent stress
*/
double plasticity_2d_vonMises_equivalent_stress(const double stress[3]);

/**
    Check if stress is inside the Von Mises yield surface.
    
    \param[in] stress stress vector [σ_x, σ_y, τ_xy]
    \param[in] sigma_y yield stress
    \return 1 if inside yield surface (elastic), 0 if outside (plastic)
*/
int plasticity_2d_vonMises_check_yield(const double stress[3], double sigma_y);

#if defined(__cplusplus)
}
#endif

#endif /* PLASTICITY_2D_VONMISES_PROJECTION_H */
