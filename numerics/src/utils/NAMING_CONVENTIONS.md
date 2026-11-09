# Numerics Naming Conventions

This document describes the standardized naming conventions applied to the Siconos numerics library.

## Overview

To improve code readability and maintainability, the following naming conventions have been established and applied across the numerics codebase.

## Standardized Terminology

### Subsystem Identifiers

| Old/Inconsistent | Standard | Description |
|-----------------|----------|-------------|
| `contact`, `cone` | `block` | Generic term for contact/cone/subsystem (preferred) |
| `contact` | `contact` | 3D friction contact (FC3D only) |
| `cone` | `cone` | Mohr-Coulomb cone (MC2D only) |
| `i`, `j`, `k` | `block_id` | Index of current subsystem |
| `nc` | `n_blocks` or `nc` | Number of subsystems (keep `nc` for backward compat) |

### Primal Variable (Reaction/Stress)

| Old/Inconsistent | Standard | Description |
|-----------------|----------|-------------|
| `reaction` | `reaction` | Global reaction vector (keep for backward compat) |
| `R` | `z` | Short form for local reaction |
| `localreaction` | `z_local` | Local subsystem reaction |
| `oldreaction` | `z_old` | Previous iteration value |

### Dual Variable (Velocity/Strain)

| Old/Inconsistent | Standard | Description |
|-----------------|----------|-------------|
| `velocity` | `velocity` | Global velocity vector (keep for backward compat) |
| - | `x` | Short form for local velocity |
| - | `x_local` | Local subsystem velocity |

### Solver Options

| Old/Inconsistent | Standard | Description |
|-----------------|----------|-------------|
| `localsolver_options` | `local_opts` | Local/internal solver options |
| `options->dparam[0]` | `SOLVER_TOL(options)` | User tolerance (macro) |
| `options->dparam[1]` | `SOLVER_RESIDUAL(options)` | Residual (macro) |
| `options->iparam[0]` | `SOLVER_MAX_ITER(options)` | Max iterations (macro) |
| `options->iparam[1]` | `SOLVER_ITER_DONE(options)` | Iterations done (macro) |

### Errors

| Convention | Description |
|------------|-------------|
| `err_full` | Full/residual error |
| `err_incr` | Incremental/light error |
| `err_rel` | Relative error |
| `err_local` | Local subsystem error |

### Iteration Counters

| Convention | Description |
|------------|-------------|
| `iter` | Current iteration number |
| `itermax` | Maximum iterations |
| `iter_done` | Iterations completed |
| `newton_iter` | Newton iteration counter |
| `ls_iter` | Line search iteration counter |

### Problem Data

| Old/Inconsistent | Standard | Description |
|-----------------|----------|-------------|
| `problem` | `problem` | Global problem structure (keep) |
| `localproblem` | `pb_local` | Local problem structure (preferred) |
| `M` | `M` | Global matrix |
| `Mlocal`, `MLocal` | `M_local` | Local subsystem matrix |
| `q` | `q` | Global right-hand side |
| `qlocal`, `qLocal` | `q_local` | Local right-hand side |

## Accessor Macros

The `naming_conventions.h` header provides standardized accessor macros:

```c
/* Get user tolerance from solver options */
SOLVER_TOL(options)

/* Get local solver tolerance */
LOCAL_SOLVER_TOL(local_opts)

/* Set local solver tolerance */
SET_LOCAL_SOLVER_TOL(local_opts, value)

/* Get residual from solver options */
SOLVER_RESIDUAL(options)

/* Set residual in solver options */
SET_SOLVER_RESIDUAL(options, value)

/* Get max iterations from solver options */
SOLVER_MAX_ITER(options)

/* Get iterations done from solver options */
SOLVER_ITER_DONE(options)

/* Set iterations done in solver options */
SET_SOLVER_ITER_DONE(options, value)

/* Get current block/contact number from solver options */
SOLVER_CURRENT_BLOCK(options)

/* Set current block/contact number in solver options */
SET_SOLVER_CURRENT_BLOCK(options, block_id)
```

## Utility Macros

```c
/* Block sizes */
BLOCK_SIZE_2D  // 2
BLOCK_SIZE_3D  // 3
BLOCK_SIZE_5D  // 5

/* Compute global index from block index and local dimension */
GLOBAL_INDEX(block_id, local_dim, local_idx)

/* Compute block start index */
BLOCK_START(block_id, local_dim)
```

## Type Aliases

```c
block_id_t    // unsigned int for block indices
iter_t        // int for iteration counters
error_t       // double for error values
tolerance_t   // double for tolerance values
```

## Files Updated

The following files have been updated with the naming conventions header:

- `numerics/src/FrictionContact/fc3d_nsgs.c`
- `numerics/src/FrictionContact/fc3d_onecontact_nonsmooth_Newton_solvers.c`
- `numerics/src/Plasticity/mc2d_nsgs.c`
- `numerics/src/Plasticity/mc2d_onecone_nonsmooth_Newton_solvers.c`
- `numerics/src/RollingFrictionContact/rolling_fc2d_nsgs.c`

## Deprecated Names

The following names are deprecated and should not be used in new code.  Backward-compatibility aliases are centralized in `numerics/src/utils/deprecated_numerics_api.h`.

### Variable naming

| Deprecated | Replacement |
|-----------|-------------|
| `contact` (as generic index) | `block_id` |
| `cone` (as generic index) | `block_id` |
| `R` (reaction) | `z` or `reaction` |
| `var_z` | `z` or `reaction` |
| `var_x` | `x` or `velocity` |
| `localsolver_options` | `local_opts` |
| `dparam[0]` | `SOLVER_TOL()` |
| `iparam[0]` | `SOLVER_MAX_ITER()` |

### Public function names

| Deprecated | Replacement |
|-----------|-------------|
| `fc3d_dsfp_set_default` | `fc3d_de_saxce_fixed_point_set_default` |
| `fc3d_fpp_set_default` | `fc3d_fixed_point_projection_set_default` |
| `fc3d_hp_set_default` | `fc3d_hyperplane_projection_set_default` |
| `fc3d_eg_set_default` | `fc3d_extra_gradient_set_default` |
| `fc3d_pfp_set_default` | `fc3d_panagiotopoulos_fixed_point_set_default` |
| `fc3d_aclmfp_set_default` | `fc3d_alart_curnier_lm_fixed_point_set_default` |
| `fc3d_tfp_set_default` | `fc3d_tresca_fixed_point_set_default` |
| `fc3d_nsn_ac_set_default` | `fc3d_nonsmooth_newton_alart_curnier_set_default` |
| `fc3d_nsn_fb_set_default` | `fc3d_nonsmooth_newton_fischer_burmeister_set_default` |
| `fc3d_nsn_nm_set_default` | `fc3d_nonsmooth_newton_jean_moreau_set_default` |
| `fc3d_nsn_ac_new_set_default` | `fc3d_nonsmooth_newton_alart_curnier_reformulated_set_default` |
| `fc3d_ipm_snm_set_default` | `fc3d_IPM_SNM_set_default` |
| `gfc3d_nsn_ac_set_default` | `gfc3d_nonsmooth_newton_alart_curnier_set_default` |
| `gfc3d_admm_set_default` | `gfc3d_alternating_direction_method_of_multipliers_set_default` |
| `gfc3d_aclmfp_set_default` | `gfc3d_alart_curnier_lm_fixed_point_set_default` |
| `gfc3d_ipm_snm_set_default` | `gfc3d_IPM_SNM_set_default` |
| `computeAlartCurnier2D_new` | `computeAlartCurnier2D_reformulated` |
| `fc3d_nonsmooth_Newton_AlartCurnier_new` | `fc3d_nonsmooth_Newton_AlartCurnier_reformulated` |
| `gfc3d_nonsmooth_Newton_AlartCurnier_new_wr` | `gfc3d_nonsmooth_Newton_AlartCurnier_reformulated_wr` |
| `NM_internalData_new` | `NM_internalData_alloc` |
| `freeAVI` | `AVI_free` |
| `freeNCP` | `NCP_free` |
| `gfc3d_vi_eg_set_default` | `gfc3d_VI_extra_gradient_set_default` |
| `gfc3d_vi_fpp_set_default` | `gfc3d_VI_fixed_point_projection_set_default` |
| `gfc2d_vi_eg_set_default` | `gfc2d_VI_extra_gradient_set_default` |
| `fc3d_FixedP_initialize` | `fc3d_fixed_point_initialize` |
| `fc3d_FixedP_solve` | `fc3d_fixed_point_solve` |
| `fc3d_FixedP_free` | `fc3d_fixed_point_free` |
| `fc3d_FixedP_computeError` | `fc3d_fixed_point_compute_error` |
| `F_GlockerFixedP` | `F_Glocker_fixed_point` |
| `compute_Z_GlockerFixedP` | `compute_Z_Glocker_fixed_point` |
| `frictionContact_display` | `frictionContactProblem_display` |
| `frictionContact_printInFile` | `frictionContactProblem_printInFile` |
| `frictionContact_printInFilename` | `frictionContactProblem_printInFilename` |
| `globalFrictionContact_display` | `globalFrictionContactProblem_display` |
| `globalFrictionContact_free` | `globalFrictionContactProblem_free` |
| `globalFrictionContact_printInFile` | `globalFrictionContactProblem_printInFile` |
| `globalFrictionContact_printInFileName` | `globalFrictionContactProblem_printInFilename` |

### Removed public functions

The following functions had no active call sites in Siconos and were removed rather than renamed:

| Removed | Rationale |
|---------|-----------|
| `Compute_NCP_error1` | Dead code; only referenced in commented-out lines |
| `Compute_NCP_error2` | Dead code; only referenced in commented-out lines |

## Future Work

- Apply expanded names to remaining cryptic abbreviations in the public API
- Remove deprecated aliases after one release cycle
- Update external examples and tutorials to use canonical names
- Apply conventions to remaining solver files
