# vkernel test suite: CppUnit integration and reference checks

Status: design approved 2026-10-28; pending spec review.
Branch: `vkernel-cppunit-tests`, forked from `main-devel-cxd-relation-sets` @ `d1ad3e4cf`.

## 1. Goal

- Every vkernel test is a CppUnit suite, registered through the project's
  test framework (`new_test`) and run by ctest.
- `bouncing_ball` gains a numeric reference check against the tutorial
  reference `BouncingBallTS.ref`, the only pre-existing external truth.
- Tests that produce a numerical trajectory gain a committed reference file
  and a regression check.

Non-goals: changing any simulated model or physics; changing ctest test
names; touching the high-level (non-vkernel) test suites.

## 2. Context

### Current inventory

| Kind | Location | Tests |
|---|---|---|
| Already CppUnit, framework-registered | `storage/tests`, `simul/tests`, `collision/tests` | `StorageTest`, `t1`, `RtOsiTest`, `rigid_shapes` |
| Legacy plain `main()`, hand-rolled `add_test` | `tests/CMakeLists.txt` | `bouncing_ball`, `bouncing_ball_registry`, `bouncing_ball_split`, `registry_remove`, `unbounded_bouncing_ball`, `bouncing_balls`, `bouncing_disk`, `bouncing_disk_line` |
| Legacy plain `main()`, hand-rolled `add_test` | `collision/tests/CMakeLists.txt` | `cxd_batch_smoke`, `polytope_cxd`, `bouncing_polytopes`, `churn_polytopes`, `polytope_fsegment`, `polytopemesh`, `polytopefdisk`, `polytopedisk`, `mixed_bodies`, `polytope_on_fem`, `bouncing_disks`, `fixed_disks`, `io_contact_points`, `sort_interactions` (registry and static_db builds) |

`fem_t3`, `bouncing_disk_sf`, `t2`, `tnonos` are built or registered only
partially today; they keep their current status except `fem_t3`, which gains
a reference check when it is registered.

### House pattern

`vkernel/soa/vkernel_tests.cmake` registers tests with

```cmake
begin_tests(src/siconos/storage/tests DEPS "vkernel;CPPUNIT::CPPUNIT")
new_test(SOURCES StorageTest.cpp ${SIMPLE_TEST_MAIN})
```

`SIMPLE_TEST_MAIN` is `kernel/tests-common/TestMain.cpp`, which runs
`CppUnit::TestFactoryRegistry::getRegistry().makeTest()`. `begin_tests`
copies `*.ref`, `*.dat`, `*.msh`, ... from the source test directory into the
build directory; `new_test` sets `TIMEOUT` and a `FAIL_REGULAR_EXPRESSION`.

### Existing verification in the legacy tests

The legacy tests are not unverified. They use a local
`check(bool, message)` helper that throws `std::runtime_error` (11 files), or
a failure counter returned as the process exit code
(`registry_remove`, `io_contact_points`). ctest therefore already fails on a
wrong result. The conversion replaces the local idiom with CppUnit
assertions; it does not add or weaken coverage.

### Reference comparison is not reusable as-is

`siconos::algebra::io::compareRefFile` takes a `SiconosDenseMatrix` and lives
in `kernel`. vkernel tests link `numerics` only, so reusing it would drag
`kernel` into every vkernel test.

### Probe: bouncing ball against the tutorial reference

`bouncing_ball` (tmax 10, h 0.005, 5 columns) against
`siconos-tutorials/.../BouncingBall/BouncingBallTS.ref` (2001 x 5, with a
leading `2001 5` header line):

- Frobenius error of the difference: `3.6e-13`
- largest element difference: `2.8e-14`
- rows differing by more than `1e-12`: none

The check is meaningful and clears the proposed tolerance by five orders of
magnitude.

## 3. Design

### 3.1 Registration

Both vkernel test directories move to `begin_tests`/`new_test` in
`vkernel/soa/vkernel_tests.cmake`. The hand-rolled `add_executable`/
`add_test` blocks in `tests/CMakeLists.txt` and
`collision/tests/CMakeLists.txt` are removed.

ctest names are preserved by passing `NAME vkernel_<test>`:

```cmake
new_test(NAME vkernel_bouncing_ball SOURCES bouncing_ball.cpp
         ${SIMPLE_TEST_MAIN} DATA bouncing_ball.ref
         DEPS "vkernel;numerics;CPPUNIT::CPPUNIT")
```

The `_sets0`/`_sets1` oracle for `bouncing_disks`, `churn_polytopes`,
`mixed_bodies`, `polytope_on_fem`, `sort_interactions` keeps its current
`add_executable` + `compare_outputs.cmake` form: it asserts that two builds
agree with each other, an invariant a single CppUnit process cannot express.
Its two executables are built from the same sources, which now contain a
CppUnit suite; the oracle compares their stdout, which stays identical.

### 3.2 Test source shape

Each legacy `main()` becomes:

```cpp
class XxxTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(XxxTest);
  CPPUNIT_TEST(testXxx);
  CPPUNIT_TEST_SUITE_END();
 public:
  void testXxx();
};
CPPUNIT_TEST_SUITE_REGISTRATION(XxxTest);
```

The model-building body moves into `testXxx()` unchanged. The local `check`
calls become `CPPUNIT_ASSERT_MESSAGE(message, condition)`. The executable
entry point comes from `SIMPLE_TEST_MAIN`; the file's own `main()` is
removed. No model code is rewritten.

### 3.3 Reference helper

New header-only, Eigen-only helper `tests/ref_check.hpp` (shared by both
test directories):

- `read_ref(path)` returns an `Eigen::MatrixXd`, skipping a leading
  non-numeric header line of the form `rows cols` (the tutorials format) and
  otherwise reading every whitespace-separated number.
- `relative_error(dat, ref)` returns
  `||dat - ref||_F / max(||ref||_F, 1.0)`, and the row and column of the
  largest element difference for the failure message.
- `CPPUNIT_ASSERT_REF(dat, "name.ref", tol)` asserts the relative error is
  below `tol`, and on failure reports both matrices' dimensions, the computed
  error, the tolerance and the worst element's position and value.
- Dimension mismatch is an explicit failure, never a partial comparison.

### 3.4 Reference data

- `bouncing_ball.ref`: copy of the tutorial `BouncingBallTS.ref`, verbatim,
  header line included. This is external truth.
- Generated references, produced by the current implementation and committed
  next to the source: `bouncing_ball_registry`, `bouncing_ball_split`,
  `unbounded_bouncing_ball`, `bouncing_balls`, `bouncing_disk`,
  `bouncing_disk_line`, `bouncing_disks`, `fixed_disks`, `fem_t3`.
  These lock current behaviour against regression; they are not independent
  truth, and the spec records that.
- Total added reference data is roughly half a megabyte.

### 3.5 Tolerance

Relative Frobenius tolerance `1e-8` for every generated reference. The
bouncing ball clears it by five orders of magnitude. A generated reference
that cannot hold `1e-8` across optimisation levels is loosened only with a
comment recording the observed spread; it is never tightened to hide noise.

### 3.6 What each test asserts

| Test | Observable | Assertion |
|---|---|---|
| `bouncing_ball` | `result.dat` 5 columns | `.ref` from the tutorial, `1e-8` |
| `bouncing_ball_registry` | `result.dat` | generated `.ref`, `1e-8` |
| `bouncing_ball_split` | `result.dat` | generated `.ref`, `1e-8` |
| `unbounded_bouncing_ball` | `result.dat` | generated `.ref`, `1e-8` |
| `bouncing_balls` | `result-many.dat` | generated `.ref`, `1e-8` |
| `bouncing_disk` | `result.dat` | generated `.ref`, `1e-8` |
| `bouncing_disk_line` | `result.dat` | generated `.ref`, `1e-8` |
| `bouncing_disks` | `result-many.dat` | generated `.ref`, `1e-8` |
| `fixed_disks` | `result-many.dat` | generated `.ref`, `1e-8` |
| `fem_t3` | `fem_t3.dat` | generated `.ref`, `1e-8` (when registered) |
| `registry_remove` | in-memory state | existing five checks as assertions: two balls remain, positions followed their items, BC set migrated, no stale keys |
| `cxd_batch_smoke` | in-memory state | existing checks as assertions |
| `polytope_cxd` | in-memory state | existing checks as assertions |
| `bouncing_polytopes` | in-memory state | existing checks as assertions |
| `churn_polytopes` | in-memory state | existing checks as assertions |
| `polytope_fsegment` | in-memory state | existing checks as assertions |
| `polytopemesh` | in-memory state | existing checks as assertions |
| `polytopedisk` | in-memory state | existing checks as assertions: contact status/h/normal, finite-difference Jacobian `1e-6`, stack heights, churn counts |
| `polytopefdisk` | in-memory state | existing checks as assertions: distances and normals, finite-difference Jacobian, resting height, no sliding or rotation |
| `mixed_bodies` | in-memory state | existing checks as assertions: relation variant per shape pair, all four categories seen, stack heights, at rest, no torque |
| `polytope_on_fem` | in-memory state | existing checks as assertions: manifold point counts, block deflection, disk/square agreement `1e-8`, mesh churn counts |
| `sort_interactions` (x2) | in-memory state, stdout | existing checks as assertions: stable sort by relation type, multiset preserved, dense keys contiguous, at least two relation types sorted |
| `io_contact_points` | in-memory state | required `y[1] == 0` and finite values as assertions, plus a non-vacuity assertion that contacts were observed |

## 4. Delivery

Two commits, each building and tested on its own:

1. **Framework and reference infrastructure, with the bouncing ball as the
   first converted test.** `ref_check.hpp`, the `new_test` registration for
   both directories, `bouncing_ball` converted with its tutorial reference,
   and the generated references for the other trajectory tests, all
   registered.
2. **The remaining structural tests converted**, one suite each, their
   `check` calls turned into assertions.

## 5. Verification

- Build in `guix shell -f guix.scm` from the worktree.
- `ctest -R '^vkernel_'` green, with the same set of test names as before the
  conversion (`ctest -N` diff must be empty except for intended additions).
- `bouncing_ball` reports a reference error near `4e-13`, and the assertion
  demonstrably fails when the reference is perturbed (checked once during
  development, not committed).
- Each generated reference is produced from a binary built at the commit that
  adds it.
- Rerun after a clean rebuild to confirm the generated references are stable.

## 6. Risks

- `new_test` sets `FAIL_REGULAR_EXPRESSION "FAILURE;Exception;ERROR;Assertion"`;
  legacy stdout containing those words would fail a passing test. Outputs are
  checked before conversion.
- `begin_tests` copies data files at configure time; adding a `.ref` requires
  a CMake re-run.
- Generated references are tied to this toolchain and machine. The tolerance
  is chosen for that, and the reference files stay revisable.
- Converting 13 collision tests changes their link line (`CPPUNIT::CPPUNIT`,
  `numerics`, `pthread`, neighbour search); each is verified by building it.
