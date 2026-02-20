# Generic NSGS (NonSmooth Gauss-Seidel) Framework

This framework provides a unified, agnostic implementation of the NSGS algorithm that can be used across different friction contact problem types (fc2d, fc3d, mc2d, rfc2d, rfc3d).

## Files Created

1. **nsgs_generic.h** - Core generic NSGS solver with callback-based architecture
2. **nsgs_generic_driver.h** - Simplified driver interface for problem-specific adapters
3. **fc3d_nsgs_generic_adapter.c** - Example adapter showing how to use with fc3d

## Architecture

The framework separates the NSGS loop logic from problem-specific details through callbacks:

### Key Structures

```c
NSGSLocalToolkit      - Contains all problem-specific callbacks
NSGSProblemData       - Common problem data (nc, q, M, mu, dimension)
NSGSProblemCallbacks  - Simplified callback structure for the driver
```

### Callback Types

- **NSGSUpdateLocalProblem** - Build/update local problem for a contact
- **NSGSSolveLocal** - Solve the local problem
- **NSGSComputeError** - Compute global error
- **NSGSLightError** - Compute light error for efficient convergence checking
- **NSGSAcceptLocal** - Accept/filter local reaction into global solution
- **NSGSAllocLocalProblem** - Allocate local problem structure
- **NSGSSetLocalTol** - Set tolerance for internal solver

## Usage

### Option 1: Using the Generic Driver (Recommended)

```c
#include "nsgs_generic_driver.h"

// Define your callbacks
NSGSProblemCallbacks callbacks = {
    .problem_name = "fc3d",
    .dimension = 3,
    .update = fc3d_update_local,
    .solve = fc3d_solve_local,
    .compute_error = fc3d_compute_error,
    .alloc_local = fc3d_alloc_local,
    .incremental_error = fc3d_incr_error,  // Optional
    .accept_local = fc3d_accept_local, // Optional
};

// Call the solver
nsgs_driver(problem, reaction, velocity, &info, options, 
            &callbacks, nc, q, M);
```

### Option 2: Using the Core Solver Directly

```c
#include "nsgs_generic.h"

// Setup problem data
NSGSProblemData problem_data = {
    .nc = problem->numberOfContacts,
    .q = problem->q,
    .M = problem->M,
    .dimension = 3
};

// Setup toolkit with callbacks
NSGSLocalToolkit toolkit = {
    .update_local_problem = my_update_function,
    .solve_local = my_solve_function,
    .compute_error = my_error_function,
    .dimension = 3,
    .use_incremental_error = 1,
    // ... other options
};

// Solve
nsgs_solve(problem, reaction, velocity, &info, options, 
           &toolkit, &problem_data);
```

## Benefits

1. **Code Reuse**: Common NSGS loop logic is written once and reused
2. **Consistency**: All problem types use the same convergence criteria, iteration counting, etc.
3. **Maintainability**: Bug fixes and improvements apply to all solvers
4. **Flexibility**: New problem types can use NSGS by implementing the callback interface
5. **Testing**: Single test suite can validate NSGS behavior across all problem types

## Migration Path

To migrate an existing solver (e.g., fc2d_nsgs):

1. Extract the local problem callbacks from the existing code
2. Wrap them with the appropriate signature
3. Create a thin wrapper that calls `nsgs_driver()`
4. Keep the original interface for backward compatibility

### Example Migration: fc2d_nsgs

```c
// Original fc2d_nsgs - becomes a thin wrapper
void fc2d_nsgs(FrictionContactProblem *problem, double *z, double *w, 
               int *info, SolverOptions *options) {
    
    // Setup callbacks
    static NSGSProblemCallbacks callbacks = {
        .problem_name = "fc2d",
        .dimension = 2,
        .update = fc2d_nsgs_update_wrapper,
        .solve = fc2d_nsgs_solve_wrapper,
        .compute_error = fc2d_compute_error_wrapper,
        .alloc_local = fc2d_alloc_local_wrapper,
        .incremental_error = fc2d_incr_error_wrapper,
    };
    
    // Call generic solver
    nsgs_driver(problem, z, w, info, options, &callbacks,
                problem->numberOfContacts, problem->q, problem->M);
}
```

## Future Enhancements

1. **Parallel NSGS**: Add OpenMP support for parallel contact processing
2. **Hybrid Solvers**: Combine NSGS with other solvers (e.g., NSGS + PATH)
3. **Adaptive Strategies**: Automatic selection of local solvers based on problem characteristics
4. **GPU Support**: CUDA/OpenCL implementations of the NSGS loop
5. **Matrix-Free**: Support for matrix-free operations using matrix-vector products

## Integration with Existing Code

The framework is designed to coexist with existing solvers. No changes to existing interfaces are required - the generic solver can be called internally by existing solver functions.
