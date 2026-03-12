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

/*!\file fc3d_short_names.h
 * \brief Short macro aliases for long FrictionContact enum names
 *
 * This header provides shorter names for frequently used enum values.
 * It is automatically included by FrictionContact headers.
 *
 * Usage: Instead of SICONOS_FRICTION_3D_NSGS, use FC3D_NSGS
 *        Instead of SICONOS_GLOBAL_FRICTION_3D_NSGS, use GFC3D_NSGS
 *
 * The short names follow the pattern:
 *   FC3D_*    = 3D Friction Contact solvers (local formulation)
 *   FC2D_*    = 2D Friction Contact solvers
 *   GFC3D_*   = Global Friction Contact 3D solvers
 *   GFC2D_*   = Global Friction Contact 2D solvers
 *   OC_*      = One-contact solvers (local solvers used within NSGS)
 */

#ifndef FC3D_SHORT_NAMES_H
#define FC3D_SHORT_NAMES_H

#include "FrictionContact_options.h"

/* ===========================================================================
 * 2D Friction Contact Solvers (FC2D_*)
 * =========================================================================== */
#define FC2D_NSGS      SICONOS_FRICTION_2D_NSGS
#define FC2D_CPG       SICONOS_FRICTION_2D_CPG
#define FC2D_LEMKE     SICONOS_FRICTION_2D_LEMKE
#define FC2D_ENUM      SICONOS_FRICTION_2D_ENUM

/* ===========================================================================
 * 3D Friction Contact Solvers - Local Formulation (FC3D_*)
 * =========================================================================== */
#define FC3D_NSGS           SICONOS_FRICTION_3D_NSGS
#define FC3D_NSGSV          SICONOS_FRICTION_3D_NSGSV
#define FC3D_PROX           SICONOS_FRICTION_3D_PROX
#define FC3D_TFP            SICONOS_FRICTION_3D_TFP
#define FC3D_NSN_AC         SICONOS_FRICTION_3D_NSN_AC
#define FC3D_NSN_AC_NEW     SICONOS_FRICTION_3D_NSN_AC_NEW
#define FC3D_NSN_FB         SICONOS_FRICTION_3D_NSN_FB
#define FC3D_NSN_NM         SICONOS_FRICTION_3D_NSN_NM
#define FC3D_DSFP           SICONOS_FRICTION_3D_DSFP
#define FC3D_VI_FPP         SICONOS_FRICTION_3D_VI_FPP
#define FC3D_VI_EG          SICONOS_FRICTION_3D_VI_EG
#define FC3D_HP             SICONOS_FRICTION_3D_HP
#define FC3D_FPP            SICONOS_FRICTION_3D_FPP
#define FC3D_EG             SICONOS_FRICTION_3D_EG
#define FC3D_GAMS_PATH      SICONOS_FRICTION_3D_GAMS_PATH
#define FC3D_GAMS_PATHVI    SICONOS_FRICTION_3D_GAMS_PATHVI
#define FC3D_GAMS_LCP_PATH  SICONOS_FRICTION_3D_GAMS_LCP_PATH
#define FC3D_GAMS_LCP_PATHVI SICONOS_FRICTION_3D_GAMS_LCP_PATHVI
#define FC3D_ACLMFP         SICONOS_FRICTION_3D_ACLMFP
#define FC3D_SOCLCP         SICONOS_FRICTION_3D_SOCLCP
#define FC3D_PFP            SICONOS_FRICTION_3D_PFP
#define FC3D_ADMM           SICONOS_FRICTION_3D_ADMM
#define FC3D_IPM            SICONOS_FRICTION_3D_IPM
#define FC3D_IPM_SNM        SICONOS_FRICTION_3D_IPM_SNM
#define FC3D_NCPG_FP        SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint
#define FC3D_NCPG_PATH      SICONOS_FRICTION_3D_NCPGlockerFBPATH
#define FC3D_NCPG_NEWTON    SICONOS_FRICTION_3D_NCPGlockerFBNewton
#define FC3D_CONVEXQP_PG_CYLINDER  SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER
#define FC3D_VI_FPP_Cylinder       SICONOS_FRICTION_3D_VI_FPP_Cylinder

/* ===========================================================================
 * One-Contact Solvers (OC_*) - Local solvers used within NSGS
 * =========================================================================== */
#define OC_NSN              SICONOS_FRICTION_3D_ONECONTACT_NSN
#define OC_NSN_GP           SICONOS_FRICTION_3D_ONECONTACT_NSN_GP
#define OC_NSN_GP_HYBRID    SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID
#define OC_PROJ             SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone
#define OC_PROJ_LI          SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration
#define OC_PROJ_REG         SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization
#define OC_PROJ_DIAG        SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization
#define OC_PROJ_V           SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity
#define OC_QUARTIC          SICONOS_FRICTION_3D_ONECONTACT_QUARTIC
#define OC_QUARTIC_NU       SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU
#define OC_CYLINDER         SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder
#define OC_CYLINDER_LI      SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration

/* ===========================================================================
 * Global Friction Contact 3D Solvers (GFC3D_*)
 * =========================================================================== */
#define GFC3D_NSGS_WR       SICONOS_GLOBAL_FRICTION_3D_NSGS_WR
#define GFC3D_NSGSV_WR      SICONOS_GLOBAL_FRICTION_3D_NSGSV_WR
#define GFC3D_PROX_WR       SICONOS_GLOBAL_FRICTION_3D_PROX_WR
#define GFC3D_DSFP_WR       SICONOS_GLOBAL_FRICTION_3D_DSFP_WR
#define GFC3D_TFP_WR        SICONOS_GLOBAL_FRICTION_3D_TFP_WR
#define GFC3D_NSGS          SICONOS_GLOBAL_FRICTION_3D_NSGS
#define GFC3D_NSN_AC_WR     SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR
#define GFC3D_NSN_AC        SICONOS_GLOBAL_FRICTION_3D_NSN_AC
#define GFC3D_GAMS_PATH     SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH
#define GFC3D_GAMS_PATHVI   SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI
#define GFC3D_VI_FPP        SICONOS_GLOBAL_FRICTION_3D_VI_FPP
#define GFC3D_VI_EG         SICONOS_GLOBAL_FRICTION_3D_VI_EG
#define GFC3D_ACLMFP        SICONOS_GLOBAL_FRICTION_3D_ACLMFP
#define GFC3D_ADMM          SICONOS_GLOBAL_FRICTION_3D_ADMM
#define GFC3D_ADMM_WR       SICONOS_GLOBAL_FRICTION_3D_ADMM_WR
#define GFC3D_IPM           SICONOS_GLOBAL_FRICTION_3D_IPM
#define GFC3D_IPM_WR        SICONOS_GLOBAL_FRICTION_3D_IPM_WR
#define GFC3D_IPM_SNM       SICONOS_GLOBAL_FRICTION_3D_IPM_SNM
#define GFC3D_IPM_SNM_WR    SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_WR
#define GFC3D_IPM_SNM_PROX  SICONOS_GLOBAL_FRICTION_3D_IPM_SNM_PROX

#endif /* FC3D_SHORT_NAMES_H */
