# Cohesive Friction Contact Tests

This directory contains tests for the cohesive friction-contact solver.

## Tests

### test_cohesive_friction_3d_simple.c

A simple test case with:
- 1 contact point
- 1 cohesive point  
- All matrices (W, V, U, X) set to identity
- Simple q_v and q_u vectors

This test verifies that the solver can:
1. Build the global M and q from block components
2. Converge to a solution satisfying the friction cone condition
3. Satisfy complementarity conditions

## Building and Running

```bash
# From build directory
cd /Users/vincent/siconos/build

# Configure with tests enabled
cmake .. -DWITH_TESTS=ON

# Build
make

# Run the test
./numerics/src/CohesiveFrictionContact/test/test_cohesive_friction_3d_simple
```

## Expected Output

The solver should converge with:
- Normal reaction r_n > 0 (contact is active)
- Tangential reaction ||r_t|| <= mu * r_n (friction cone satisfied)
- r_n * u_n ≈ 0 (complementarity)
