# Tools Directory (Legacy)

This directory now serves as a backward compatibility layer.

Most files have been moved to specialized directories:
- `core/` - Fundamental types and options (SolverOptions, etc.)
- `linalg/` - Linear algebra (NumericsMatrix, etc.)
- `optimization/` - Optimization methods (line search, etc.)
- `io/` - Input/Output (HDF5, etc.)
- `utils/` - General utilities

The headers in this directory are wrappers that include the actual headers
from their new locations. They are provided for backward compatibility only
and may be deprecated in a future release.
