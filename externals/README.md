
# Siconos externals

This directory contains sources for software libraries required or
used optionally by Siconos.  Some sources may be modified for use
within Siconos.  Copyright information is documented in a file named
"copyright" found in each first-level subfolder.

Most of these libraries are linked together into a library called
`libsiconos_externals`, which is used principally by numerics.  Not all
are: `cxd` and `boost_hana_pfr` are header-only, while `CompactNSearch`
and `cuNSearch` are built from their own `CMakeLists.txt`.  Each entry
below names the part of Siconos that uses it.

Directories `tools` and `renderer` are native to the Siconos
project.

Directories `PATH_SDK`, `sort`, and `optim_misc` contain sources that
are not permissive as defined by https://opensource.org/osd-annotated

These directories may be removed and Siconos can still compile and
function.  Currently the LCP_QP solver will give runtime errors if
`optim_misc/ql0001/ql0001.f` is removed.

Directory lbl contains the C wrapper C wrapper around MA27 and MA57
written by D. Orban. LBL is a unified C interface to the multifrontal
symmetric indefinite linear system solvers MA27 and MA57 from the Harwell
Subroutine Library (HSL). The HSL library is not permissive. To use it,
you must ask the library to HSL and put is in lbl/ext. For some compatiblity
reason, we also inlude metis 4.
The source are taken form github https://github.com/dpo/lbl
commit 80d468c1584e475287bcd3a2762a4ee3db5d2c15

Directory CompactNSearch contains the neighbour search of
https://github.com/InteractiveComputerGraphics/CompactNSearch, used by
the vkernel component in builds without GPU.
commit a9ab7c7 (`git describe` 1.0.0-42-ga9ab7c7)
Local modifications are the "[vkernel] CompactNSearch: ..." commits,
listed by `git log -- externals/CompactNSearch`.

Directory cuNSearch contains the CUDA neighbour search of
https://github.com/InteractiveComputerGraphics/cuNSearch, used by the
vkernel component when WITH_GPU is ON.
commit 4dd5598
Local modifications are the "[vkernel] cuNSearch: ..." commits, listed
by `git log -- externals/cuNSearch`.

Directory boost_hana_pfr contains the header-only part of Boost that the
vkernel component uses (hana and pfr).  It is placed under its own name
rather than under `boost`, so that it cannot be mistaken for a complete
Boost tree; the directory itself is on the include path, so that
`#include "boost/hana.hpp"` and `#include "boost/pfr.hpp"` keep the
upstream spelling.  The Boost version it was taken from is not recorded.
It is unmodified, and it must come before the system Boost on the
include path when nvcc is used (see the comment in vkernel/CMakeLists.txt).

Directory cxd contains CXD, a convex collision detection library (GJK,
EPA, and SAT for small meshes), used by the vkernel component only.  No
upstream URL is recorded: the sources are maintained here, under the
Apache License 2.0 (see cxd/COPYING).  In practice it is header-only, the
`cxd` target of cxd/src/CMakeLists.txt not being added to the build.
