
# Siconos externals

This directory contains sources for software libraries required or
used optionally by Siconos.  Some sources may be modified for use
within Siconos.  Copyright information is documented in a file named
"copyright" found in each first-level subfolder.

These libraries are linked together into a library called
`libsiconos_externals`, which is used principally by numerics, but may
also be needed in other parts of Siconos.

Directories `tools`, `swig`, and `renderer` are native to the Siconos
project.

Directories `PATH_SDK`, `sort`, and `optim_misc` contain sources that
are not permissive as defined by https://opensource.org/osd-annotated

These directories may be removed and Siconos can still compile and
function.  Currently the LCP_QP solver will give runtime errors if
`optim_misc/ql0001/ql0001.f` is removed.
