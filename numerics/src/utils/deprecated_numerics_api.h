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

/*!\file deprecated_numerics_api.h
 * \brief Centralized backward-compatibility aliases for renamed numerics functions.
 *
 * This header is intentionally separate from the main API headers so that the
 * canonical code base is not polluted with old names.  Downstream code that
 * still relies on the old names can include this header once and continue to
 * compile while migrating to the new names.
 *
 * The aliases below are deprecated and may be removed in a future release.
 */

#ifndef DEPRECATED_NUMERICS_API_H
#define DEPRECATED_NUMERICS_API_H

#include "NumericsFwd.h"

/* Portable deprecation attribute. */
#if defined(__GNUC__) || defined(__clang__)
#define SICONOS_DEPRECATED(msg) __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
#define SICONOS_DEPRECATED(msg) __declspec(deprecated(msg))
#else
#define SICONOS_DEPRECATED(msg)
#endif

/* ============================================================================
 * FrictionContact 3D solver set_default aliases
 * ============================================================================ */

SICONOS_DEPRECATED("Use fc3d_de_saxce_fixed_point_set_default instead")
static inline void fc3d_dsfp_set_default(SolverOptions* options) {
  fc3d_de_saxce_fixed_point_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_fixed_point_projection_set_default instead")
static inline void fc3d_fpp_set_default(SolverOptions* options) {
  fc3d_fixed_point_projection_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_hyperplane_projection_set_default instead")
static inline void fc3d_hp_set_default(SolverOptions* options) {
  fc3d_hyperplane_projection_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_extra_gradient_set_default instead")
static inline void fc3d_eg_set_default(SolverOptions* options) {
  fc3d_extra_gradient_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_panagiotopoulos_fixed_point_set_default instead")
static inline void fc3d_pfp_set_default(SolverOptions* options) {
  fc3d_panagiotopoulos_fixed_point_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_alart_curnier_lm_fixed_point_set_default instead")
static inline void fc3d_aclmfp_set_default(SolverOptions* options) {
  fc3d_alart_curnier_lm_fixed_point_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_tresca_fixed_point_set_default instead")
static inline void fc3d_tfp_set_default(SolverOptions* options) {
  fc3d_tresca_fixed_point_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_nonsmooth_newton_alart_curnier_set_default instead")
static inline void fc3d_nsn_ac_set_default(SolverOptions* options) {
  fc3d_nonsmooth_newton_alart_curnier_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_nonsmooth_newton_fischer_burmeister_set_default instead")
static inline void fc3d_nsn_fb_set_default(SolverOptions* options) {
  fc3d_nonsmooth_newton_fischer_burmeister_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_nonsmooth_newton_jean_moreau_set_default instead")
static inline void fc3d_nsn_nm_set_default(SolverOptions* options) {
  fc3d_nonsmooth_newton_jean_moreau_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_nonsmooth_newton_alart_curnier_reformulated_set_default instead")
static inline void fc3d_nsn_ac_new_set_default(SolverOptions* options) {
  fc3d_nonsmooth_newton_alart_curnier_reformulated_set_default(options);
}

SICONOS_DEPRECATED("Use fc3d_IPM_SNM_set_default instead")
static inline void fc3d_ipm_snm_set_default(SolverOptions* options) {
  fc3d_IPM_SNM_set_default(options);
}

/* ============================================================================
 * Global FrictionContact 3D solver set_default aliases
 * ============================================================================ */

SICONOS_DEPRECATED("Use gfc3d_nonsmooth_newton_alart_curnier_set_default instead")
static inline void gfc3d_nsn_ac_set_default(SolverOptions* options) {
  gfc3d_nonsmooth_newton_alart_curnier_set_default(options);
}

SICONOS_DEPRECATED("Use gfc3d_alternating_direction_method_of_multipliers_set_default instead")
static inline void gfc3d_admm_set_default(SolverOptions* options) {
  gfc3d_alternating_direction_method_of_multipliers_set_default(options);
}

SICONOS_DEPRECATED("Use gfc3d_alart_curnier_lm_fixed_point_set_default instead")
static inline void gfc3d_aclmfp_set_default(SolverOptions* options) {
  gfc3d_alart_curnier_lm_fixed_point_set_default(options);
}

SICONOS_DEPRECATED("Use gfc3d_IPM_SNM_set_default instead")
static inline void gfc3d_ipm_snm_set_default(SolverOptions* options) {
  gfc3d_IPM_SNM_set_default(options);
}

/* ============================================================================
 * Misleading _new suffix aliases
 * ============================================================================ */

SICONOS_DEPRECATED("Use computeAlartCurnier2D_reformulated instead")
static inline void computeAlartCurnier2D_new(double* reaction, double* velocity, double mu,
                                             double* rho, double* residual, double* jacobian_u,
                                             double* jacobian_r) {
  computeAlartCurnier2D_reformulated(reaction, velocity, mu, rho, residual, jacobian_u,
                                     jacobian_r);
}

SICONOS_DEPRECATED("Use fc3d_nonsmooth_Newton_AlartCurnier_reformulated instead")
static inline void fc3d_nonsmooth_Newton_AlartCurnier_new(FrictionContactProblem* problem,
                                                          double* reaction, double* velocity,
                                                          int* info, SolverOptions* options) {
  fc3d_nonsmooth_Newton_AlartCurnier_reformulated(problem, reaction, velocity, info, options);
}

SICONOS_DEPRECATED("Use gfc3d_nonsmooth_Newton_AlartCurnier_reformulated_wr instead")
static inline void gfc3d_nonsmooth_Newton_AlartCurnier_new_wr(
    GlobalFrictionContactProblem* problem, double* reaction, double* velocity,
    double* globalVelocity, int* info, SolverOptions* options) {
  gfc3d_nonsmooth_Newton_AlartCurnier_reformulated_wr(problem, reaction, velocity,
                                                      globalVelocity, info, options);
}

SICONOS_DEPRECATED("Use NM_internalData_alloc instead")
static inline void NM_internalData_new(NumericsMatrix* A) { NM_internalData_alloc(A); }

/* ============================================================================
 * Problem struct lifecycle verb aliases
 * ============================================================================ */

SICONOS_DEPRECATED("Use AVI_free instead")
static inline void freeAVI(AffineVariationalInequalities* avi) { AVI_free(avi); }

SICONOS_DEPRECATED("Use NCP_free instead")
static inline void freeNCP(NonlinearComplementarityProblem* ncp) { NCP_free(ncp); }

/* ============================================================================
 * Wave 2: VI solver set_default aliases
 * ============================================================================ */

SICONOS_DEPRECATED("Use gfc3d_VI_extra_gradient_set_default instead")
static inline void gfc3d_vi_eg_set_default(SolverOptions* options) {
  gfc3d_VI_extra_gradient_set_default(options);
}

SICONOS_DEPRECATED("Use gfc3d_VI_fixed_point_projection_set_default instead")
static inline void gfc3d_vi_fpp_set_default(SolverOptions* options) {
  gfc3d_VI_fixed_point_projection_set_default(options);
}

SICONOS_DEPRECATED("Use gfc2d_VI_extra_gradient_set_default instead")
static inline void gfc2d_vi_eg_set_default(SolverOptions* options) {
  gfc2d_VI_extra_gradient_set_default(options);
}

/* ============================================================================
 * Wave 2: Glocker fixed-point local solver aliases
 * ============================================================================ */

SICONOS_DEPRECATED("Use fc3d_fixed_point_initialize instead")
static inline void fc3d_FixedP_initialize(FrictionContactProblem* problem,
                                          FrictionContactProblem* localproblem,
                                          SolverOptions* localsolver_option) {
  fc3d_fixed_point_initialize(problem, localproblem, localsolver_option);
}

SICONOS_DEPRECATED("Use fc3d_fixed_point_solve instead")
static inline int fc3d_FixedP_solve(FrictionContactProblem* localproblem, double* reaction,
                                    SolverOptions* options) {
  return fc3d_fixed_point_solve(localproblem, reaction, options);
}

SICONOS_DEPRECATED("Use fc3d_fixed_point_free instead")
static inline void fc3d_FixedP_free(FrictionContactProblem* localproblem,
                                    SolverOptions* localsolver_option) {
  fc3d_fixed_point_free(localproblem, localsolver_option);
}

SICONOS_DEPRECATED("Use fc3d_fixed_point_compute_error instead")
static inline double fc3d_FixedP_computeError(int contact, int dimReaction, double* reaction,
                                              double* error) {
  return fc3d_fixed_point_compute_error(contact, dimReaction, reaction, error);
}

SICONOS_DEPRECATED("Use F_Glocker_fixed_point instead")
static inline void F_GlockerFixedP(int sizeF, double* reaction, double* FVector, int up2Date) {
  F_Glocker_fixed_point(sizeF, reaction, FVector, up2Date);
}

SICONOS_DEPRECATED("Use compute_Z_Glocker_fixed_point instead")
static inline void compute_Z_GlockerFixedP(int contact, double* reactionstep) {
  compute_Z_Glocker_fixed_point(contact, reactionstep);
}

/* ============================================================================
 * Wave 2: FrictionContactProblem lifecycle verb aliases
 * ============================================================================ */

SICONOS_DEPRECATED("Use frictionContactProblem_display instead")
static inline void frictionContact_display(FrictionContactProblem* problem) {
  frictionContactProblem_display(problem);
}

SICONOS_DEPRECATED("Use frictionContactProblem_printInFile instead")
static inline int frictionContact_printInFile(FrictionContactProblem* problem, FILE* file) {
  return frictionContactProblem_printInFile(problem, file);
}

SICONOS_DEPRECATED("Use frictionContactProblem_printInFilename instead")
static inline int frictionContact_printInFilename(FrictionContactProblem* problem,
                                                  char* filename) {
  return frictionContactProblem_printInFilename(problem, filename);
}

/* ============================================================================
 * Wave 2: GlobalFrictionContactProblem lifecycle verb aliases
 * ============================================================================ */

SICONOS_DEPRECATED("Use globalFrictionContactProblem_display instead")
static inline void globalFrictionContact_display(GlobalFrictionContactProblem* problem) {
  globalFrictionContactProblem_display(problem);
}

SICONOS_DEPRECATED("Use globalFrictionContactProblem_free instead")
static inline void globalFrictionContact_free(GlobalFrictionContactProblem* problem) {
  globalFrictionContactProblem_free(problem);
}

SICONOS_DEPRECATED("Use globalFrictionContactProblem_printInFile instead")
static inline int globalFrictionContact_printInFile(GlobalFrictionContactProblem* problem,
                                                    FILE* file) {
  return globalFrictionContactProblem_printInFile(problem, file);
}

SICONOS_DEPRECATED("Use globalFrictionContactProblem_printInFilename instead")
static inline int globalFrictionContact_printInFileName(GlobalFrictionContactProblem* problem,
                                                        char* filename) {
  return globalFrictionContactProblem_printInFilename(problem, filename);
}

#undef SICONOS_DEPRECATED

#endif /* DEPRECATED_NUMERICS_API_H */
