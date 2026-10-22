/* Scale benchmark for the cxd-backed polytope narrowphase.

   SCOPE: NARROWPHASE ONLY. No time-stepper, no contact solver, no
   broadphase. The solver stack cannot run at scale: one_step_integrator's
   assemble_setup() resizes w_matrix_assembled to a DENSE
   raw_inter_size x raw_inter_size matrix (raw_inter_size = 2 x number of
   interactions), so 40k interactions would ask for an 80k x 80k matrix
   (51 GB). mass_matrix_assembled is dense too (3N x 3N). This benchmark
   therefore drives the narrowphase stage on its own, which is the only
   part of the collision pipeline that scales today.

   Everything runs at step column 0. With no one_step_integrator in the
   world there is no keep<> on q/velocity/y/ydot, so every step index maps
   to column 0 and ds.q() is the current geometry.

   What is measured, per step:
     - the whole stage (rebind + checksum + pack + batch + write-back),
     - the same split phase by phase (the stage's own timers),
     - the analytic disk-disk reference on the SAME pair set, as the
       yardstick the polytope path has to be compared against. That
       reference is a plain arithmetic loop: it carries NO storage-handle
       overhead, so the ratio is an upper bound on what the batch path
       would cost against a hand-written analytic polytope path.

   Usage: polytope_narrowphase_scale [N [steps [spacing jitter]]]
          without N the benchmark sweeps N = 250, 1000, 4000. */

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <print>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "siconos/siconos.hpp"

// clang-format off
/* the relation + stage headers must be visible before the config names
   them; the load-bearing cxd include order lives inside cxd_narrowphase.hpp */
#include "siconos/collision/polytopepolytope_r.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;

struct polytope : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using polytope_shape = collision::shape::polytope;
using diskdisk_r = collision::diskdisk_r;

/* topo is the acyclic marker of the relation/narrowphase pair; defined
   right after the interaction that names it */
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;

struct interaction : simul::interaction<nslaw, diskdisk_r, polytopepolytope_r> {};
struct topo : simul::topology<polytope, interaction> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

/* No simulation and no one_step_integrator: the world holds only what the
   narrowphase stage touches. */
struct data_t
    : storage::make<
          env, topo, polytope, polytope_shape, diskdisk_r, nslaw, interaction,
          polytopepolytope_r, cxd_narrowphase,
          storage::with_properties<
              storage::wrapped<polytope, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<interaction, some::unbounded_collection>,
              storage::wrapped<polytopepolytope_r, some::unbounded_collection>,
              storage::attached<polytope, storage::pattern::symbol<"shape">,
                                storage::some::item_ref<polytope_shape>>,
              storage::time_invariant<storage::pattern::attr_t<polytope, "fext">>,
              storage::diagonal<storage::pattern::attr_t<polytope, "mass_matrix">>>> {};
}  // namespace siconos::config

namespace {

using namespace siconos;
using namespace storage;

/* assert() is compiled out under Release: an explicit throwing check keeps
   the RED/GREEN honest */
void check(bool ok, const std::string& what) {
  if (!ok) throw std::runtime_error("CHECK FAILED: " + what);
}

/* xorshift64*: a fixed-seed PRNG, so the benchmark world is reproducible
   without depending on the implementation-defined std:: engines */
struct rng {
  std::uint64_t s;
  explicit rng(std::uint64_t seed) : s(seed ? seed : 0x9E3779B97F4A7C15ull) {}
  std::uint64_t next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * 2685821657736338717ull;
  }
  double uniform(double a, double b) {
    return a + (b - a) * (static_cast<double>(next() >> 11) / 9007199254740992.0);
  }
};

/* the body shape: an axis-aligned square of side kSize, centred on the
   body origin (the same square the bouncing_polytopes test uses) */
constexpr double kSize = 2.0;
constexpr double kCircumradius = 1.4142135623730951;  // sqrt(2)

struct layout {
  std::size_t n = 0;
  std::size_t side = 0;  // grid side
  double pitch = 0.0;    // grid pitch
  std::vector<double> x0, y0, vx, vy;
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
};

/* N bodies on a ceil(sqrt N) squared grid, jittered, with four neighbours
   per cell (east, north, north-east, north-west) so that every pair is
   emitted exactly once, from its lower cell */
layout build_layout(std::size_t n, double spacing_factor, double jitter_factor) {
  layout L;
  L.n = n;
  L.side = static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(n))));
  L.pitch = spacing_factor * kSize;
  rng r(0xC0FFEEull);
  L.x0.resize(n);
  L.y0.resize(n);
  L.vx.resize(n);
  L.vy.resize(n);
  for (std::size_t gy = 0; gy < L.side; ++gy) {
    for (std::size_t gx = 0; gx < L.side; ++gx) {
      const std::size_t i = gy * L.side + gx;
      if (i >= n) break;
      L.x0[i] =
          static_cast<double>(gx) * L.pitch + r.uniform(-1.0, 1.0) * jitter_factor * kSize;
      L.y0[i] =
          static_cast<double>(gy) * L.pitch + r.uniform(-1.0, 1.0) * jitter_factor * kSize;
      L.vx[i] = r.uniform(-1.0, 1.0) * 0.01;
      L.vy[i] = r.uniform(-1.0, 1.0) * 0.01;
    }
  }
  for (std::size_t gy = 0; gy < L.side; ++gy) {
    for (std::size_t gx = 0; gx < L.side; ++gx) {
      const std::size_t i = gy * L.side + gx;
      if (i >= n) break;
      constexpr int dx[4] = {+1, 0, +1, -1};
      constexpr int dy[4] = {0, +1, +1, +1};
      for (int k = 0; k < 4; ++k) {
        const long nx = static_cast<long>(gx) + dx[k];
        const long ny = static_cast<long>(gy) + dy[k];
        if (nx < 0 || ny < 0) continue;
        if (nx >= static_cast<long>(L.side) || ny >= static_cast<long>(L.side)) continue;
        const std::size_t j =
            static_cast<std::size_t>(ny) * L.side + static_cast<std::size_t>(nx);
        if (j >= n) continue;
        L.pairs.emplace_back(std::min(i, j), std::max(i, j));
      }
    }
  }
  return L;
}

/* the shape shared by every body */
auto make_shape(config::data_t& data) {
  using namespace siconos;
  using namespace storage;

  auto shape = storage::add<config::polytope_shape>(data);
  const std::array<std::array<double, 2>, 4> square = {{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
  for (int i = 0; i < 4; ++i) shape.vertices()[i] = {square[i][0], square[i][1]};
  shape.nvertices() = 4;
  return shape;
}

/* The narrowphase must reproduce the analytic gap of two axis-aligned
   squares: |y2 - y1| - kSize. This is the reference the bouncing_polytopes
   test pins, run here on its own world so that the benchmark fails loudly
   instead of printing meaningless timings. */
void self_check() {
  using namespace siconos;
  using namespace storage;

  config::data_t data;

  auto shape = make_shape(data);
  auto topo_h = storage::add<config::topo>(data);
  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = 0.9;
  auto dsds = storage::add<config::diskdisk_r>(data);

  auto b1 = storage::add<config::polytope>(data);
  auto b2 = storage::add<config::polytope>(data);
  for (auto b : {b1, b2}) {
    b.q() = {0.0, 0.0, 0.0};
    b.velocity() = {0.0, 0.0, 0.0};
    b.mass_matrix().diagonal() << 1.0, 1.0, 2.0 / 3.0;
    b.fext() = {0.0, 0.0, 0.0};
    storage::prop<"shape">(b) = shape;
  }
  b2.q() = {0.0, 3.0, 0.0};

  auto inter = topo_h.link(b1, b2);
  inter.nslaw() = nslaw;
  inter.relation() = dsds;

  auto np = storage::add<config::cxd_narrowphase>(data);
  np.dsds_r() = dsds;
  np.ensure(0);

  auto rels = storage::handles<config::polytopepolytope_r>(data, 0);
  std::size_t n_relations = 0;
  double h = 0.0;
  std::size_t status = 0;
  for (auto rel : rels) {
    ++n_relations;
    h = rel.h();
    status = rel.status();
  }
  check(n_relations == 1, "the stage rebound exactly one relation");
  check(std::isfinite(h), "the cached gap is finite");
  check(std::abs(h - 1.0) <= 1e-9, "the cached gap is the analytic gap (3 - 2 = 1)");
  check(status == static_cast<std::size_t>(cxd::status::separated),
        "the pair is reported as separated");
}

struct measurement {
  std::size_t n = 0;
  std::size_t pairs = 0;
  std::size_t contacts = 0;
  std::size_t enclosed = 0;
  double ms_total = 0.0;  // per step
  double ms_rebind = 0.0;
  double ms_checksum = 0.0;
  double ms_pack = 0.0;
  double ms_batch = 0.0;
  double ms_writeback = 0.0;
  double ns_per_pair = 0.0;
  double ns_per_pair_disk = 0.0;
  double ratio = 0.0;
};

measurement run_scale(std::size_t n, std::size_t n_steps, double spacing_factor,
                      double jitter_factor) {
  using namespace siconos;
  using namespace storage;

  config::data_t data;
  const layout L = build_layout(n, spacing_factor, jitter_factor);

  auto shape = make_shape(data);
  auto topo_h = storage::add<config::topo>(data);
  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = 0.9;
  auto dsds = storage::add<config::diskdisk_r>(data);

  /* The bodies are added ONE BY ONE and the handles are kept. The count
     overload of storage::add returns a LAZY view whose transform calls
     add() again on every walk, so a range walked inside the step loop would
     grow the world at every step; the items are materialised once here. */
  using body_handle_t = decltype(storage::add<config::polytope>(data));
  std::vector<body_handle_t> bodies;
  bodies.reserve(L.n);
  for (std::size_t i = 0; i < L.n; ++i) bodies.push_back(storage::add<config::polytope>(data));

  for (std::size_t i = 0; i < bodies.size(); ++i) {
    auto& body = bodies[i];
    body.q() = {L.x0[i], L.y0[i], 0.0};
    body.velocity() = {L.vx[i], L.vy[i], 0.0};
    body.mass_matrix().diagonal() << 1.0, 1.0, 2.0 / 3.0;
    body.fext() = {0.0, 0.0, 0.0};
    storage::prop<"shape">(body) = shape;
  }

  for (const auto& [i, j] : L.pairs) {
    auto inter = topo_h.link(
        storage::make_handle(data, storage::index<config::polytope, std::size_t>(i)),
        storage::make_handle(data, storage::index<config::polytope, std::size_t>(j)));
    inter.nslaw() = nslaw;
    inter.relation() = dsds;
  }

  auto np = storage::add<config::cxd_narrowphase>(data);
  np.dsds_r() = dsds;

  measurement m;
  m.n = n;
  m.pairs = L.pairs.size();

  /* the first ensure() creates the relations and the batch state */
  np.ensure(0);
  check(np.n_pairs() == m.pairs, "the stage packed every pair");

  const double dt = 0.005;
  std::vector<double> xf(L.x0), yf(L.y0);
  double seconds_ensure = 0.0;
  for (std::size_t s = 0; s < n_steps; ++s) {
    /* move every body: the world checksum must change, otherwise ensure()
       short-circuits on the cache and the batch never re-runs */
    const double a = dt * static_cast<double>(s + 1);
    for (std::size_t i = 0; i < bodies.size(); ++i) {
      xf[i] = L.x0[i] + L.vx[i] * a;
      yf[i] = L.y0[i] + L.vy[i] * a;
      bodies[i].q() = {xf[i], yf[i], 0.0};
    }
    const auto t0 = std::chrono::steady_clock::now();
    np.ensure(0);
    const auto t1 = std::chrono::steady_clock::now();
    seconds_ensure += std::chrono::duration<double>(t1 - t0).count();
    m.ms_rebind += np.t_rebind() * 1e3;
    m.ms_checksum += np.t_checksum() * 1e3;
    m.ms_pack += np.t_pack() * 1e3;
    m.ms_batch += np.t_batch() * 1e3;
    m.ms_writeback += np.t_writeback() * 1e3;
  }

  m.contacts = np.n_contacts();
  m.enclosed = np.n_enclosed();
  const double nsteps = static_cast<double>(n_steps);
  m.ms_total = seconds_ensure * 1e3 / nsteps;
  m.ms_rebind /= nsteps;
  m.ms_checksum /= nsteps;
  m.ms_pack /= nsteps;
  m.ms_batch /= nsteps;
  m.ms_writeback /= nsteps;

  /* the analytic disk-disk reference on the SAME pair set, at the SAME
     final positions: a plain arithmetic loop, no storage handles */
  double sink = 0.0;
  const auto d0 = std::chrono::steady_clock::now();
  for (const auto& [i, j] : L.pairs) {
    sink += std::hypot(xf[j] - xf[i], yf[j] - yf[i]) - 2.0 * kCircumradius;
  }
  const auto d1 = std::chrono::steady_clock::now();
  check(std::isfinite(sink), "the analytic reference is finite");

  m.ns_per_pair = m.ms_total * 1e6 / static_cast<double>(m.pairs);
  m.ns_per_pair_disk =
      std::chrono::duration<double>(d1 - d0).count() * 1e9 / static_cast<double>(m.pairs);
  m.ratio = m.ns_per_pair / m.ns_per_pair_disk;
  return m;
}

}  // namespace

int main(int argc, char* argv[]) {
  using namespace siconos;

  std::size_t n_bodies = 0;  // 0 -> sweep
  std::size_t n_steps = 20;
  double spacing_factor = 1.2;
  double jitter_factor = 0.3;

  if (argc > 1) n_bodies = std::strtoul(argv[1], nullptr, 10);
  if (argc > 2) n_steps = std::strtoul(argv[2], nullptr, 10);
  if (argc > 3) spacing_factor = std::atof(argv[3]);
  if (argc > 4) jitter_factor = std::atof(argv[4]);
  if (n_steps == 0) n_steps = 1;

  try {
    self_check();
  } catch (const std::exception& e) {
    std::print("self check FAILED: {}\n", e.what());
    return 1;
  }

  std::vector<std::size_t> sizes;
  if (n_bodies)
    sizes.push_back(n_bodies);
  else
    sizes = {250, 1000, 4000};

  std::print("polytope narrowphase scale benchmark\n");
  std::print(
      "  scope: narrowphase only (no time-stepper, no contact solver, no broadphase)\n");
  std::print(
      "  world: convex polytopes (square, side {}), 4-neighbour grid, spacing {} x size,"
      " jitter {} x size, {} steps\n",
      kSize, spacing_factor, jitter_factor, n_steps);
  std::print(
      "  the disk column is a plain arithmetic loop on the same pair set: no"
      " storage-handle overhead\n");
  std::print(
      "\n{:>7} {:>8} {:>9} {:>7} {:>10} {:>9} {:>9} {:>8} {:>9} {:>9} {:>10} {:>13} {:>7}\n",
      "N", "pairs", "contacts", "encl", "ms/step", "rebind", "checksum", "pack", "batch",
      "writeback", "ns/pair", "disk ns/pr", "ratio");

  for (const auto n : sizes) {
    const auto m = run_scale(n, n_steps, spacing_factor, jitter_factor);
    std::print(
        "{:>7} {:>8} {:>9} {:>7} {:>10.3} {:>9.3} {:>9.3} {:>8.3} {:>9.3} {:>9.3} {:>10.1} "
        "{:>13.1} {:>7.1}\n",
        m.n, m.pairs, m.contacts, m.enclosed, m.ms_total, m.ms_rebind, m.ms_checksum,
        m.ms_pack, m.ms_batch, m.ms_writeback, m.ns_per_pair, m.ns_per_pair_disk, m.ratio);
  }
  return 0;
}
