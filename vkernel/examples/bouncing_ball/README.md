# bouncing_ball: an external project against an installed vkernel

This is not part of the in-tree build and it is not installed. It exists to
prove one thing: after `cmake --install`, a program that has never seen the
Siconos source tree compiles and runs against the installed vkernel, through
the include strings

```c++
#include "siconos/vkernel.hpp"
#include "siconos/vkernel/config/config_builder.hpp"
```

and the CMake target `Siconos::vkernel`.

It is the same model as `vkernel/src/siconos/vkernel/tests/bouncing_ball.cpp`
with the test harness removed: a ball of mass 1 kg, radius 0.1 m, dropped from
1 m, bouncing for 10 s at h = 0.005 with restitution 0.9 under the Moreau-Jean
scheme.

## Build and run

Install Siconos first, then, with `<prefix>` the directory it was configured
with (in this worktree `<workdir>/install-clang`):

```bash
cmake -S vkernel/examples/bouncing_ball -B /tmp/vk-example \
      -DCMAKE_PREFIX_PATH=<prefix>
cmake --build /tmp/vk-example -j4
LD_LIBRARY_PATH=<prefix>/lib /tmp/vk-example/bouncing_ball
```

Expected: exit 0 and one line

```
t = 10.000000  q = ...  v = ...
```

with a finite, non-zero `q`.

`LD_LIBRARY_PATH` is needed because this executable is not installed and has
no RPATH.

## Two things that do not work, and why

- `find_package(Siconos COMPONENTS vkernel)` fails: nothing sets
  `siconos_vkernel_FOUND`, so `check_required_components` rejects the request.
  Use plain `find_package(Siconos REQUIRED)`.
- `cmake --install . --prefix <elsewhere>` breaks the package: the generated
  `siconos-config.cmake` bakes the configure-time prefix into
  `set_and_check(siconos_ROOT ...)`, so a package installed elsewhere points at
  a directory that does not exist. Install into the configured prefix.
