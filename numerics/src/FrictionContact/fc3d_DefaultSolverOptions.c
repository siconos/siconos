/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdio.h>                                     // for size_t
#include "Friction_cst.h"                              // for SICONOS_FRICTI...
#include "Newton_methods.h"                            // for newton_lsa_set...
#include "NumericsFwd.h"                               // for SolverOptions
#include "SolverOptions.h"                             // for SolverOptions
#include "fc3d_Solvers.h"                              // 
#include "fc2d_Solvers.h"                              // 
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"        // for fc3d_nonsmooth...
#include "fc3d_nonsmooth_Newton_FischerBurmeister.h"   // for fc3d_nonsmooth...
#include "fc3d_nonsmooth_Newton_natural_map.h"         // for fc3d_nonsmooth...
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"  // for fc3d_onecontac...
#include "fc3d_projection.h"                           // for fc3d_projectio...
#include "fc3d_unitary_enumerative.h"                  // for fc3d_unitary_e...

int fc3d_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  solver_options_nullify(options);

  int info = -1;
  switch (solverId)
    {

      // FC3D Solvers
    case SICONOS_FRICTION_3D_NSGS:
      {
        info =    fc3d_nsgs_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_NSGSV:
      {
        info =    fc3d_nsgs_velocity_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ADMM:
      {
        info =    fc3d_admm_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_PROX:
      {
        info =    fc3d_proximal_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_TFP:
      {
        info =    fc3d_TrescaFixedPoint_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_PFP:
      {
        info =    fc3d_Panagiotopoulos_FixedPoint_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ACLMFP:
      {
        info =    fc3d_ACLMFixedPoint_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_SOCLCP:
      {
        info =    fc3d_SOCLCP_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_DSFP:
      {
        info =    fc3d_DeSaxceFixedPoint_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_FPP:
      {
        info =    fc3d_fixedPointProjection_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_EG:
      {
        info =    fc3d_ExtraGradient_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_VI_FPP:
      {
        info =    fc3d_VI_FixedPointProjection_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_VI_EG:
      {
        info =    fc3d_VI_ExtraGradient_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_HP:
      {
        info =    fc3d_HyperplaneProjection_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_NSN_AC:
      {
        info =    fc3d_nonsmooth_Newton_AlartCurnier_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_NSN_AC_TEST:
      {
        size_t iSize = 20;
        size_t dSize = 20;
        size_t iter_max = 200;
        double tol = 1e-12;
        solver_options_fill(options, solverId, iSize, dSize, iter_max, tol);
        newton_lsa_setDefaultSolverOptions(options);
        options->solverId = solverId;
        //    options->iparam[5] = 1;
        //    options->iparam[7] = 1;
        options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] = SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED;
        /* 0 STD AlartCurnier, 1 JeanMoreau, 2 STD generated, 3 JeanMoreau generated */

        break;
      }
    case SICONOS_FRICTION_3D_NSN_FB:
      {
        info =    fc3d_nonsmooth_Newton_FischerBurmeister_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_NSN_NM:
      {
        info =    fc3d_nonsmooth_Newton_NaturalMap_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder:
      {
        info =    fc3d_ConvexQP_ProjectedGradient_Cylinder_setDefaultSolverOptions(options);
        options->solverId = SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder;
        break;
      }
    case SICONOS_FRICTION_3D_VI_FPP_Cylinder:
      {
        info =    fc3d_VI_FixedPointProjection_Cylinder_setDefaultSolverOptions(options);
        options->solverId = SICONOS_FRICTION_3D_VI_FPP_Cylinder;
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC:
      {
        info =    fc3d_unitary_enumerative_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU:
      {
        info =    fc3d_unitary_enumerative_setDefaultSolverOptions(options);
        options->solverId = SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU;
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization:
      {
        info =    fc3d_projectionOnConeWithDiagonalization_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization:
      {
        info =    fc3d_projectionOnConeWithRegularization_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone:
      {
        info =    fc3d_projectionOnCone_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
      {
        info =    fc3d_projectionOnConeWithLocalIteration_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity:
      {
        info =    fc3d_projectionOnCone_velocity_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder:
      {
        info =    fc3d_projectionOnCylinder_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration:
      {
        info =    fc3d_projectionOnCylinderWithLocalIteration_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_NSN:
      {
        info =  fc3d_onecontact_nonsmooth_Newton_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP:
      {
        info =  fc3d_onecontact_nonsmooth_Newton_gp_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID:
      {
        info =  fc3d_onecontact_nonsmooth_Newton_gp_setDefaultSolverOptions(options);
        options->solverId=SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID;
        break;
      }

      // FC2D solvers
    case SICONOS_FRICTION_2D_NSGS:
      {
        info =    fc2d_sparse_nsgs_setDefaultSolverOptions(options);
        /* info =    fc2d_nsgs_setDefaultSolverOptions(options); */
        break;
      }
    case SICONOS_FRICTION_2D_PGS:
      {
        info =    fc2d_nsgs_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_2D_CPG:
      {
        info =    fc2d_cpg_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_2D_LATIN:
      {
        info =    fc2d_latin_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_2D_LEMKE:
      {
        info =    fc2d_lexicolemke_setDefaultSolverOptions(options);
        break;
      }
    case SICONOS_FRICTION_2D_ENUM:
      {
        info =    fc2d_enum_setDefaultSolverOptions(options);
        break;
      }

    default:
      {
        solver_options_set(options, solverId);

      }
    }

  return info;
}
