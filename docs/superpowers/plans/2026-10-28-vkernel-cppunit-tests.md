# vkernel tests: CppUnit and reference checks — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development
> (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run every vkernel test as a CppUnit suite registered through `new_test`, with a
numeric reference check on every test that produces a trajectory.

**Architecture:** Two commits. Commit (a) adds the header-only `ref_check.hpp`, converts the
seven trajectory tests in `src/siconos/tests/` and the two in `collision/tests/` to CppUnit
suites, registers them through `begin_tests`/`new_test` (ctest names unchanged), and commits a
reference file for each. Commit (b) converts the remaining structural tests and turns their
local `check()` calls into CppUnit assertions.

**Tech Stack:** C++23, clang++, Eigen, CppUnit, CMake (`new_test`), ctest, guix.

**Spec:** `docs/superpowers/specs/2026-10-28-vkernel-cppunit-tests-design.md`

## Global Constraints

- Commit messages: one line, at most 68 characters.
- Exactly two commits for the delivery, each building and tested on its own.
- Build inside `guix shell -f guix.scm -- <cmd>`; never outside it.
- Run test binaries with `LD_LIBRARY_PATH=/tmp/opencode/extra-libs`.
- Never touch the `main-devel-cxd-relation-sets` worktree.
- ctest names stay `vkernel_<test>`; `ctest -N -R '^vkernel_'` must list the same names.
- Reference tolerance is relative Frobenius `1e-8`; loosen only with a comment recording the
  observed spread, never tighten.
- Style: `.clang-format`, ColumnLimit 95; pre-commit runs clang-format on `*.cpp/*.hpp`
  (`.dat` and `.ref` are excluded).
- C++ standard for the tests is `--std=c++23` (`std::print`); `new_test` does not apply
  `vkernel/soa`'s `${flags}`, so each converted target gets
  `set_target_properties(<t> PROPERTIES COMPILE_FLAGS "${flags}")`, as `RtOsiTest` and `t1`
  already do.

## Spec deltas (decided while planning, recorded here and in the spec)

1. **`fem_t3` cannot be registered.** It reads `square_200.msh`, which is not in the repo, and
   aborts (`Cannot open file: square_200.msh`). It keeps its current status: built, not
   registered, no `.ref`. Its row in the spec's assertion table stays a future item.
2. **`bouncing_disks` writes no file.** `std::fstream result_file("result-many.dat")` opens
   without creating the file, so in a fresh working directory nothing is written. One token
   changes to `std::ofstream`; no model change.
3. **`bouncing_balls` writes one row.** The per-step `std::print` in its loop goes to stdout,
   not to `result_file`, so `result-many.dat` holds only the initial row. One token changes to
   `std::print(result_file, ...)`; no model change.
4. **Helper location.** `vkernel/soa/src/siconos/tests/ref_check.hpp`, included as
   `"ref_check.hpp"` from `tests/` and `"../../tests/ref_check.hpp"` from `collision/tests/`
   (quote include is relative to the including file, so no CMake include path is added).

## Review Focus

- A `.ref` that the test cannot find: the path is relative to the working directory, which is
  `build/vkernel/soa/src/siconos/<dir>` for `new_test`; a typo must fail loudly, not skip.
- A dimension mismatch between `result.dat` and the committed `.ref`: explicit failure, never a
  partial comparison.
- A test whose stdout contains `FAILURE;Exception;ERROR;Assertion` (the `new_test`
  `FAIL_REGULAR_EXPRESSION`): would fail a passing test. Checked on all 25 legacy binaries
  before conversion; only `registry_remove` prints `0 failure(s)`, which does not match
  (case-sensitive) and disappears when it becomes a suite.
- A reference that only holds at one optimization level: rerun after a clean rebuild.
- The `_sets0`/`_sets1` oracle binaries read the same `.ref` from their own working directory;
  it must be copied there, or they fail for the wrong reason.

---

### Task 1: Baseline build and records

**Files:** none (build tree only: `/tmp/siconos-build/vkernel-cppunit-tests`,
install prefix `/tmp/siconos-install/vkernel-cppunit-tests`).

**Interfaces:**
- Produces: `/tmp/siconos-build/vkernel-cppunit-tests/in_guix.sh` (runs a command in
  `guix shell -f <worktree>/guix.scm`), `configure.sh`, `build.sh`, `test.sh`, and the baseline
  files `/tmp/opencode/vkernel-ctest-names-before.txt`, `/tmp/opencode/vkernel-ctest-before.log`.

- [ ] **Step 1: Write the three scripts**

`in_guix.sh`: `exec guix shell -f <worktree>/guix.scm -- bash -c "$*"`.
`configure.sh`: `in_guix.sh "CXX=clang++ cmake <src> -DWITH_GUIX=1 -Dpetsc_ROOT=${GUIX_ENVIRONMENT}
-DHDF5_ROOT=${GUIX_ENVIRONMENT}/lib -DCOMPONENTS='externals;numerics;kernel;control;mechanics;io;vkernel'
-DCMAKE_BUILD_TYPE=Release -DWITH_PYTHON_WRAPPER=0 -DWITH_PYB11_WRAPPER=1 -DWITH_TESTING=1
-DWITH_kernel_TESTING=0 -DWITH_io_TESTING=0 -DWITH_vkernel_TESTING=1 -DWITH_numerics_TESTING=0
-DSICONOS_CUSTOM_INSTALL=/tmp/siconos-install/vkernel-cppunit-tests -DISOLATED_INSTALL=True
-DWITH_FCLIB=1 -DFCLIB_ROOT=1 -DWITH_MPI=0 -DWITH_GPU=0 -DWITH_CUDA=0 -DWITH_BULLET=1
-DWITH_HDF5=1 -DWITH_PETSC=1 -DWITH_OPENMP=1 -DCMAKE_EXPORT_COMPILE_COMMANDS=ON"`.

- [ ] **Step 2: Configure and build**

Run: `./configure.sh && ./build.sh` (build target `vkernel-tests` is not enough; build `all`).
Expected: configure and build succeed.

- [ ] **Step 3: Record the baseline**

Run: `cd /tmp/siconos-build/vkernel-cppunit-tests && ctest -N -R '^vkernel_' | sed 's/.*: //' | sort > /tmp/opencode/vkernel-ctest-names-before.txt`
Expected: 29 names (`vkernel_` prefix, 24 tests + 5 `_sets_oracle`).

Run: `LD_LIBRARY_PATH=/tmp/opencode/extra-libs ctest -R '^vkernel_' --output-on-failure | tee /tmp/opencode/vkernel-ctest-before.log`
Expected: 29/29 pass.

- [ ] **Step 4: Verify the perturbations the spec requires (not committed)**

Run the tutorial reference against the current `bouncing_ball` output
(`examples/mechanics/BouncingBall/BouncingBallTS.ref` of siconos-tutorials):
expected Frobenius error ≈ `3.6e-13`.

---

### Task 2: `ref_check.hpp`

**Files:**
- Create: `vkernel/soa/src/siconos/tests/ref_check.hpp`

**Interfaces:**
- Produces:
  - `namespace siconos::test`:
    - `Eigen::MatrixXd read_ref(const std::string& path)` — all whitespace-separated numbers
      of a file as a matrix; a leading line of exactly two integer tokens (`rows cols`) is
      skipped; the number of columns is taken from the first data line.
    - `struct RefError { double relative; Eigen::Index row; Eigen::Index col; double worst; }`
      and `RefError relative_error(const Eigen::MatrixXd& dat, const Eigen::MatrixXd& ref)` —
      `relative = ||dat - ref||_F / max(||ref||_F, 1.0)`, `worst`/`row`/`col` the largest
      absolute element difference and where it is.
  - Macro `CPPUNIT_ASSERT_REF(dat, ref_file, tol)`, usable inside any test method: asserts the
    dimensions match (explicit failure, no partial comparison), then `relative <= tol`; the
    failure message carries both dimensions, the computed error, the tolerance and the worst
    element with its 1-based row and column.

- [ ] **Step 1: Write the header**

Include `<cppunit/TestAssert.h>`, `<Eigen/Dense>`, `<fstream>`, `<sstream>`, `<string>`,
`<vector>`. `read_ref` reads all tokens into `std::vector<double>`, then reshapes with the
column count of its first line; an empty file or a ragged row throws
`std::runtime_error("<path>: ...")`. The macro evaluates `dat` once into a local reference,
calls `read_ref(ref_file)` once, and builds both messages with `std::ostringstream`. Dimension
mismatch uses `CPPUNIT_ASSERT_MESSAGE` before any arithmetic.

- [ ] **Step 2: Check it compiles and bites**

Not a ctest test yet: compile it through the first converted test (Task 3), then during Task 6
perturb one generated reference by `1e-3` and confirm the test fails with the dimension and the
worst element reported, then restore the file.

---

### Task 3: Trajectory tests in `tests/` as CppUnit suites

**Files (modify):**
- `vkernel/soa/src/siconos/tests/bouncing_ball.cpp`
- `vkernel/soa/src/siconos/tests/bouncing_ball_registry.cpp`
- `vkernel/soa/src/siconos/tests/bouncing_ball_split.cpp`
- `vkernel/soa/src/siconos/tests/unbounded_bouncing_ball.cpp`
- `vkernel/soa/src/siconos/tests/bouncing_balls.cpp`
- `vkernel/soa/src/siconos/tests/bouncing_disk.cpp`
- `vkernel/soa/src/siconos/tests/bouncing_disk_line.cpp`

**Interfaces:**
- Consumes: `siconos::test::read_ref`, `CPPUNIT_ASSERT_REF` (Task 2).
- Produces: one registered suite per file, class names `BouncingBallTest`,
  `BouncingBallRegistryTest`, `BouncingBallSplitTest`, `UnboundedBouncingBallTest`,
  `BouncingBallsTest`, `BouncingDiskTest`, `BouncingDiskLineTest`; each with a single method
  `void test<Name>()` holding the former `main()` body.

- [ ] **Step 1: Convert each file the same way**

Add `#include <cppunit/TestAssert.h>`, `#include <cppunit/extensions/HelperMacros.h>`,
`#include "ref_check.hpp"`. Replace `int main(int argc, char* argv[]) {` with

```cpp
namespace {
class XxxTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(XxxTest);
  CPPUNIT_TEST(testXxx);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testXxx() {
```

and the closing `}` with `}\n};\n\nCPPUNIT_TEST_SUITE_REGISTRATION(XxxTest);\n\n}  // namespace`.
Delete `return` statements of the former `main`. Nothing else in the body changes, except:

- `bouncing_balls.cpp:191`: `std::print(` → `std::print(result_file, ` (spec delta 3).
- each file: `result_file.close();` before the reference check, then
  `CPPUNIT_ASSERT_REF(siconos::test::read_ref("<file>"), "<name>.ref", 1e-8);` where `<file>` is
  the file the test writes (`result.dat`, or `result-many.dat` for `bouncing_balls`) and
  `<name>` is the test name.

- [ ] **Step 2: Commit nothing yet; the CMake registration is Task 5**

---

### Task 4: `bouncing_disks` and `fixed_disks` in `collision/tests/`

**Files (modify):**
- `vkernel/soa/src/siconos/collision/tests/bouncing_disks.cpp`
- `vkernel/soa/src/siconos/collision/tests/fixed_disks.cpp`
- `vkernel/soa/src/siconos/collision/tests/CMakeLists.txt` (oracle wiring only)

**Interfaces:**
- Consumes: `siconos::test::read_ref`, `CPPUNIT_ASSERT_REF` (Task 2).
- Produces: suites `BouncingDisksTest`, `FixedDisksTest`, each registered once; the same class
  serves both the ctest binary and the `_sets0`/`_sets1` oracle binaries.

- [ ] **Step 1: Convert both files as in Task 3 Step 1**, with
  `#include "../../tests/ref_check.hpp"`, `std::fstream` → `std::ofstream` at
  `bouncing_disks.cpp:189` (spec delta 2), and a `1e-8` check of `result-many.dat` against
  `bouncing_disks.ref` / `fixed_disks.ref`.

- [ ] **Step 2: Give the oracle binaries a main and the reference**

In the `foreach(test ...)` oracle block: add `${SIMPLE_TEST_MAIN}` to the `add_executable`
sources of `${test}_sets${sets}`, and copy the test's reference into the oracle working
directory at configure time so both binaries find it:

```cmake
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test}_sets_oracle)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test}.ref)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${test}.ref
    ${CMAKE_CURRENT_BINARY_DIR}/${test}_sets_oracle/${test}.ref COPYONLY)
endif()
```

Both binaries print the CppUnit runner line (`OK (1 tests)`) identically, so the stdout
comparison in `compare_outputs.cmake` is unaffected.

---

### Task 5: Framework registration for commit (a)

**Files (modify):**
- `vkernel/soa/vkernel_tests.cmake`
- `vkernel/soa/src/siconos/tests/CMakeLists.txt`
- `vkernel/soa/src/siconos/collision/tests/CMakeLists.txt`
- `docs/superpowers/specs/2026-10-28-vkernel-cppunit-tests-design.md` (spec deltas)

**Interfaces:**
- Consumes: the suite names of Tasks 3 and 4.
- Produces: ctest names `vkernel_bouncing_ball`, `vkernel_bouncing_ball_registry`,
  `vkernel_bouncing_ball_split`, `vkernel_unbounded_bouncing_ball`, `vkernel_bouncing_balls`,
  `vkernel_bouncing_disk`, `vkernel_bouncing_disk_line`, `vkernel_bouncing_disks`,
  `vkernel_fixed_disks`, unchanged from today; `SIMPLE_TEST_MAIN` is the entry point.

- [ ] **Step 1: Register them in `vkernel_tests.cmake`**

After the existing `storage/tests` block, add

```cmake
  if(WITH_GPU)
    set(COMPACTNSEARCH cuNSearch)
  else()
    set(COMPACTNSEARCH CompactNSearch)
  endif()

  begin_tests(src/siconos/tests DEPS
              "vkernel;numerics;CPPUNIT::CPPUNIT;Eigen3::Eigen;pthread")
  new_test(NAME vkernel_bouncing_ball SOURCES bouncing_ball.cpp ${SIMPLE_TEST_MAIN}
           DATA bouncing_ball.ref)
  ... one new_test per trajectory test ...
```

each with `DATA <name>.ref`; `bouncing_disk` and `bouncing_disk_line` add
`DEPS "vkernel;numerics;CPPUNIT::CPPUNIT;Eigen3::Eigen;pthread;${COMPACTNSEARCH}"`. Then a
second `begin_tests(src/siconos/collision/tests DEPS "...;${COMPACTNSEARCH}")` with
`new_test(NAME vkernel_bouncing_disks ...)` and `new_test(NAME vkernel_fixed_disks ...)`.
After the registrations, apply the compiler flags used by the old targets:

```cmake
  foreach(t vkernel_bouncing_ball ... vkernel_fixed_disks)
    set_target_properties(${t} PROPERTIES COMPILE_FLAGS "${flags}")
  endforeach()
```

- [ ] **Step 2: Delete the hand-rolled blocks for those tests**

`tests/CMakeLists.txt`: remove the eight `add_executable`, their `set_target_properties`,
`target_link_libraries` and the `add_test` foreach; keep `fem_t3` as an executable target and
the comment explaining why it is not registered.
`collision/tests/CMakeLists.txt`: drop `bouncing_disks` and `fixed_disks` from the
`add_executable`/`set_target_properties`/`target_link_libraries`/`add_test` lists; keep the
oracle block.

- [ ] **Step 3: Generate the committed references**

Configure and build, then for each of the nine generated references run the freshly built
binary in a scratch directory (`/tmp/opencode/gen/<name>`), copy its `.dat` to
`vkernel/soa/src/siconos/<dir>/<name>.ref`, and re-run configure so `begin_tests`/`DATA` copies
it into the build tree. `bouncing_ball.ref` is instead the verbatim copy of
`examples/mechanics/BouncingBall/BouncingBallTS.ref` of siconos-tutorials
(header line included).

- [ ] **Step 4: Record the spec deltas** in §3.4/§3.6 of the design document (fem_t3 not
  registered; the two one-token recording fixes) and add a note that the plan is committed
  with them.

---

### Task 6: Verify commit (a)

- [ ] **Step 1:** `./build.sh` → succeeds.
- [ ] **Step 2:** `ctest -N -R '^vkernel_' | sed 's/.*: //' | sort | diff /tmp/opencode/vkernel-ctest-names-before.txt -`
  → empty.
- [ ] **Step 3:** `LD_LIBRARY_PATH=/tmp/opencode/extra-libs ctest -R '^vkernel_' --output-on-failure`
  → 29/29 pass; `vkernel_bouncing_ball` prints a reference error near `4e-13` on failure paths.
- [ ] **Step 4:** Perturb `bouncing_ball.ref` by `1e-3` in the build tree only (not the source),
  run the test, confirm it fails and reports dimensions/error/tolerance/worst element.
- [ ] **Step 5:** Clean rebuild of the test targets, rerun the suite → still green (reference
  stability).
- [ ] **Step 6:** Commit (a).

---

### Task 7: Structural tests (commit b)

**Files (modify):** `registry_remove.cpp`, `cxd_batch_smoke.cpp`, `polytope_cxd.cpp`,
`bouncing_polytopes.cpp`, `churn_polytopes.cpp`, `polytope_fsegment.cpp`, `polytopemesh.cpp`,
`polytopedisk.cpp`, `polytopefdisk.cpp`, `mixed_bodies.cpp`, `polytope_on_fem.cpp`,
`io_contact_points.cpp`, `sort_interactions.cpp`, `rigid_shapes.cpp` (drop its own `main()`),
both test `CMakeLists.txt` (drop the remaining hand-rolled blocks) and `vkernel_tests.cmake`
(register the remaining 14 names, including `vkernel_sort_interactions_registry` /
`_static_db` with `COMPILE_DEFINITIONS SORT_STATIC_DB` for the latter and `-UNDEBUG` plus
`ENVIRONMENT OMP_NUM_THREADS=1` on both).

**Interfaces:**
- Consumes: the conversion shape of Tasks 3/4.
- Produces: suites `RegistryRemoveTest`, `CxdBatchSmokeTest`, `PolytopeCxdTest`,
  `BouncingPolytopesTest`, `ChurnPolytopesTest`, `PolytopeFsegmentTest`, `PolytopeMeshTest`,
  `PolytopeDiskTest`, `PolytopeFDiskTest`, `MixedBodiesTest`, `PolytopeOnFemTest`,
  `IoContactPointsTest`, `SortInteractionsTest`.

- [ ] **Step 1:** Wrap each `main()` in a suite (Task 3 Step 1), replacing the local
  `check(condition, message)` with `CPPUNIT_ASSERT_MESSAGE(message, condition)`; in
  `registry_remove` and `io_contact_points` turn the failure counter into assertions (and add
  `io_contact_points`' non-vacuity assertion that contacts were observed).
- [ ] **Step 2:** Register the rest in `vkernel_tests.cmake`, delete the hand-rolled blocks,
  keep the five `_sets_oracle` blocks.
- [ ] **Step 3:** Rebuild, `ctest -N` diff unchanged, full `ctest -R '^vkernel_'` green.
- [ ] **Step 4:** Commit (b), one line, ≤68 characters.

## Self-review notes

- Spec §3.6 is covered test by test; `fem_t3` is the one row not implemented, for the reason in
  the deltas.
- Every step names the exact file, name and value; no step asks the implementer to choose a
  tolerance, a path or a ctest name.
- Type names are fixed once and reused: `<Name>Test` / `test<Name>()`, `siconos::test::read_ref`,
  `CPPUNIT_ASSERT_REF`.
