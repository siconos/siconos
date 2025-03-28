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

/** \file SiconosNumerics_Solvers.h
\brief Define Macros used for solver name <---> id conversion.
See functions in SolverOptions.h
 */

// List of solvers generated with
//   git grep '^ *SICONOS_[A-Z_]* = [0-9]' | sed
//   's/.*\(SICONOS_[A-Z_]*\).*/SICONOS_SOLVER_MACRO(\1); \\/' | grep -v NUMERICS_PROBLEM >
//   NonSmoothSolvers/SiconosNumerics_Solvers.h

#undef SICONOS_SOLVER_MACRO
#define SICONOS_REGISTER_SOLVERS()                                                          \
  SICONOS_SOLVER_MACRO(SICONOS_AVI_CAOFERRIS);                                              \
  SICONOS_SOLVER_MACRO(SICONOS_AVI_PATHAVI);                                                \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_LEMKE);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_NSGS_SBM);                                               \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_PGS);                                                    \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_CPG);                                                    \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_LATIN);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_LATIN_W);                                                \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_QP);                                                     \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_NSQP);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_NEWTONMIN);                                              \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_NEWTON_FB_FBLSA);                                        \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_PSOR);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_RPGS);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_PATH);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_ENUM);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_AVI_CAOFERRIS);                                          \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_PIVOT);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_BARD);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_MURTY);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_NEWTON_MIN_FBLSA);                                       \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_PATHSEARCH);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_PIVOT_LUMOD);                                            \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_GAMS);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_LCP_CONVEXQP_PG);                                            \
  SICONOS_SOLVER_MACRO(SICONOS_MCP_OLD_FB);                                                 \
  SICONOS_SOLVER_MACRO(SICONOS_MCP_NEWTON_FB_FBLSA);                                        \
  SICONOS_SOLVER_MACRO(SICONOS_MCP_NEWTON_MIN_FBLSA);                                       \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_PGS);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_RPGS);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_PSOR);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_RPSOR);                                                 \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_PATH);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_ENUM);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_SIMPLEX);                                               \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_ENUM);                                           \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_PATH_ENUM);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_SIMPLEX);                                        \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_PATH);                                           \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_PATH_ENUM);                                      \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_FB);                                                    \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_DIRECT_FB);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_PGS_SBM);                                               \
  SICONOS_SOLVER_MACRO(SICONOS_MLCP_LCP_LEMKE);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_NCP_NEWTON_FB_FBLSA);                                        \
  SICONOS_SOLVER_MACRO(SICONOS_NCP_NEWTON_MIN_FBLSA);                                       \
  SICONOS_SOLVER_MACRO(SICONOS_NCP_PATHSEARCH);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_NCP_PATH);                                                   \
  SICONOS_SOLVER_MACRO(SICONOS_RELAY_PGS);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_RELAY_ENUM);                                                 \
  SICONOS_SOLVER_MACRO(SICONOS_RELAY_PATH);                                                 \
  SICONOS_SOLVER_MACRO(SICONOS_RELAY_LEMKE);                                                \
  SICONOS_SOLVER_MACRO(SICONOS_RELAY_AVI_CAOFERRIS);                                        \
  SICONOS_SOLVER_MACRO(SICONOS_RELAY_AVI_CAOFERRIS_TEST);                                   \
  SICONOS_SOLVER_MACRO(SICONOS_VI_EG);                                                      \
  SICONOS_SOLVER_MACRO(SICONOS_VI_FPP);                                                     \
  SICONOS_SOLVER_MACRO(SICONOS_VI_HP);                                                      \
  SICONOS_SOLVER_MACRO(SICONOS_VI_BOX_QI);                                                  \
  SICONOS_SOLVER_MACRO(SICONOS_VI_BOX_AVI_LSA);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_2D_NSGS);                                           \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_2D_CPG);                                            \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_2D_ENUM);                                           \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_2D_LEMKE);                                          \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NSGS);                                           \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NSGSV);                                          \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_PROX);                                           \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_TFP);                                            \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_PFP);                                            \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NSN_AC);                                         \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NSN_AC_TEST);                                    \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NSN_FB);                                         \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NSN_NM);                                         \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_DSFP);                                           \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_VI_FPP);                                         \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_VI_EG);                                          \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_HP);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint);                         \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_FPP);                                            \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_EG);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_NSN);                                 \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);                              \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID);                       \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization); \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone);                    \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);  \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization);  \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder);                \
  SICONOS_SOLVER_MACRO(                                                                     \
      SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration);               \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity);           \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NCPGlockerFBNewton);                             \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_NCPGlockerFBPATH);                               \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER);                           \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_GAMS_PATH);                                      \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_GAMS_PATHVI);                                    \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_GAMS_LCP_PATH);                                  \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_GAMS_LCP_PATHVI);                                \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_SOCLCP);                                         \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ACLMFP);                                         \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);                             \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU);                          \
  SICONOS_SOLVER_MACRO(SICONOS_FRICTION_3D_ADMM);                                           \
  SICONOS_SOLVER_MACRO(SICONOS_ROLLING_FRICTION_3D_NSGS);                                   \
  SICONOS_SOLVER_MACRO(                                                                     \
      SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration);           \
  SICONOS_SOLVER_MACRO(SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone);            \
  SICONOS_SOLVER_MACRO(SICONOS_ROLLING_FRICTION_3D_ADMM);                                   \
  SICONOS_SOLVER_MACRO(SICONOS_ROLLING_FRICTION_2D_NSGS);                                   \
  SICONOS_SOLVER_MACRO(                                                                     \
      SICONOS_ROLLING_FRICTION_2D_ONECONTACT_ProjectionOnConeWithLocalIteration);           \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);                                 \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR);                                \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_PROX_WR);                                 \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_DSFP_WR);                                 \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_TFP_WR);                                  \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_NSGS);                                    \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR);                               \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_NSN_AC);                                  \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH);                               \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI);                             \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_VI_EG);                                   \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_VI_FPP);                                  \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_ACLMFP);                                  \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_ADMM);                                    \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_ADMM_WR);                                 \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_FRICTION_3D_IPM);                                     \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR);                         \
  SICONOS_SOLVER_MACRO(SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM);                             \
  SICONOS_SOLVER_MACRO(SICONOS_SOCLCP_NSGS);                                                \
  SICONOS_SOLVER_MACRO(SICONOS_SOCLCP_VI_FPP);                                              \
  SICONOS_SOLVER_MACRO(SICONOS_SOCLCP_VI_EG);                                               \
  SICONOS_SOLVER_MACRO(SICONOS_CONVEXQP_VI_EG);                                             \
  SICONOS_SOLVER_MACRO(SICONOS_CONVEXQP_VI_FPP);                                            \
  SICONOS_SOLVER_MACRO(SICONOS_CONVEXQP_PG);                                                \
  SICONOS_SOLVER_MACRO(SICONOS_CONVEXQP_ADMM);                                              \
  SICONOS_SOLVER_MACRO(SICONOS_GENERIC_MECHANICAL_NSGS);                                    \
  SICONOS_SOLVER_MACRO(SICONOS_NEWTON_LSA);
