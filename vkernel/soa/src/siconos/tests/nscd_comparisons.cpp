/* The comparisons of the NSCD scheme (Moreau-Jean) with the Anitescu-Tasora
   (projected Gauss-Seidel and APGD) and explicit DEM schemes, on the cases
   of the design note ("Plan: where NSCD wins"), with the scenes of
   collision/tests/nscd_world.hpp. One CSV file per case in the output
   directory; tests/nscd_comparisons_plot.py makes the tables and the plots.

     stability      a column of 5 disks at rest, h from 1e-4 to 1e-2: the
                    largest kinetic energy and overlap (DEM is stable for h
                    below ~ 2 sqrt(m / k_n))
     column         a column of N disks, at rest and landing from gaps of
                    0.01: the overlaps (DEM: the bottom one ~ N m g / k_n)
     incline_stick  a square on an incline below the friction angle: x(t)
     sliding        a square launched at v0 with friction: the stopping
                    distance and the lift-off vs h (Anitescu-Tasora lifts
                    off, u_N = mu |u_T|)
     incline_slide  a square on an incline above the friction angle: x(T)
                    vs g (sin alpha - mu cos alpha) T^2 / 2, vs h
     sliding_path, incline_slide_path
                    the trajectories of the square in these two cases, at
                    h = 5e-3 (DEM 1e-4)

   The cost of each run: its wall time and the iterations of its solver
   (summed over the steps; none for DEM). Moreau-Jean and Anitescu-Tasora
   at the same solver tolerance, 1e-12: with 1e-8, the inexact velocities
   let the gap of a resting contact drift above the activation threshold
   (1e-10) in a few hundred steps, the contact drops out, the body falls by
   g h^2 and impacts (columns of 20 disks with Moreau-Jean, of 1 disk with
   APGD): an artifact of the tolerance, not of a scheme. DEM with k_n = 1e5, a damping ratio
   of 0.5 on a contact with the ground, k_t = 2/7 k_n, theta = 1.

   usage: nscd_comparisons <output directory> [case...] (all by default) */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <print>
#include <string>
#include <vector>

#include "siconos/collision/tests/nscd_world.hpp"

namespace {

using namespace nscd;

constexpr double kn = 1e5;
constexpr double kt = 2. / 7. * kn;
const double cn = 2 * 0.5 * std::sqrt(m * kn);
constexpr double tolerance = 1e-12;

struct scheme_t {
  const char* name;
  std::size_t id;
};
const std::vector<scheme_t> all_schemes = {{"Moreau-Jean", schemes::local_problem},
                                           {"AT (PGS)", schemes::cone_complementarity},
                                           {"AT (APGD)", schemes::cone_complementarity_apgd},
                                           {"DEM", schemes::explicit_penalty}};

bool is_dem(const scheme_t& s) { return s.id == schemes::explicit_penalty; }

scene scene_for(const scheme_t& s, double mu, double h, double tmax, double alpha = 0.) {
  scene sc{.mu = mu, .h = h, .tmax = tmax, .alpha = alpha, .scheme = s.id};
  sc.tolerance = tolerance;
  if (is_dem(s)) {
    sc.theta = 1.;
    sc.kn = kn;
    sc.cn = cn;
    sc.kt = kt;
  }
  return sc;
}

// the cost of a run, while it steps
struct cost {
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
  std::size_t steps = 0;
  std::size_t iterations = 0;
  double wall() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  }
};

// one step of the world, counted
void step(world& w, cost& c) {
  w.step();
  ++c.steps;
  auto osnspb = w.simulation.one_step_nonsmooth_problem();
  c.iterations += siconos::storage::attr<"number_of_iterations">(osnspb);
}

// a column of n disks on the ground, gap between them and above the ground
struct column_result {
  cost c;
  bool finite = true;
  double max_kinetic_energy = 0;
  double max_overlap = 0;  // the largest overlap of any contact, over the run
  double bottom_overlap_end = 0;
};

column_result column(const scheme_t& s, std::size_t n, double h, double tmax, double gap) {
  using namespace siconos;
  world w(scene_for(s, 0.3, h, tmax));
  std::vector<decltype(w.add_disk(0., 0.))> disks;
  for (std::size_t i = 0; i < n; ++i)
    disks.push_back(w.add_disk(0., r + gap + (2 * r + gap) * i));
  w.finalize();
  column_result res;
  while (w.simulation.has_next_event()) {
    step(w, res.c);
    const auto st = w.simulation.current_step();
    double below = 0., energy = 0.;
    for (std::size_t i = 0; i < n; ++i) {
      const double y = storage::attr<"q">(disks[i], st)(1);
      const auto& v = storage::attr<"velocity">(disks[i], st);
      const double overlap = (below + r) - y;
      if (i == 0) res.bottom_overlap_end = overlap;
      res.max_overlap = std::max(res.max_overlap, overlap);
      energy += 0.5 * m * (v(0) * v(0) + v(1) * v(1)) + 0.25 * m * r * r * v(2) * v(2);
      below = y + r;
    }
    if (!std::isfinite(energy) || energy > 1e6) {
      // blown up: no need to go on
      res.finite = false;
      res.max_kinetic_energy = INFINITY;
      break;
    }
    res.max_kinetic_energy = std::max(res.max_kinetic_energy, energy);
  }
  return res;
}

std::ofstream open_csv(const std::filesystem::path& dir, const std::string& name,
                       const std::string& header) {
  std::ofstream f(dir / (name + ".csv"));
  f << header << "\n";
  return f;
}

void stability(const std::filesystem::path& dir) {
  auto f = open_csv(dir, "stability",
                    "scheme,h,steps,wall_s,iterations,finite,max_kinetic_energy,max_overlap,"
                    "bottom_overlap_end");
  for (const auto& s : all_schemes) {
    for (double h : {1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2}) {
      const auto res = column(s, 5, h, 1., 0.);
      std::println(f, "{},{},{},{},{},{},{},{},{}", s.name, h, res.c.steps, res.c.wall(),
                   res.c.iterations, res.finite, res.max_kinetic_energy, res.max_overlap,
                   res.bottom_overlap_end);
      std::println(stderr, "stability {} h {:.0e}: {}", s.name, h,
                   res.finite ? "stable" : "blown up");
    }
  }
}

void columns(const std::filesystem::path& dir) {
  auto f = open_csv(dir, "column",
                    "start,scheme,n,h,steps,wall_s,iterations,max_overlap,bottom_overlap_end,"
                    "n_m_g_over_kn");
  for (const auto& [start, gap] : {std::pair{"rest", 0.}, std::pair{"landing", 0.01}}) {
    for (const auto& s : all_schemes) {
      // DEM at a tenth of the period of a contact, the others at 5e-3
      const double h = is_dem(s) ? 1e-4 : 5e-3;
      for (std::size_t n : {1, 2, 5, 10, 20, 50}) {
        const auto res = column(s, n, h, 1., gap);
        std::println(f, "{},{},{},{},{},{},{},{},{},{}", start, s.name, n, h, res.c.steps,
                     res.c.wall(), res.c.iterations, res.max_overlap, res.bottom_overlap_end,
                     n * m * g / kn);
        std::println(stderr, "column {} {} n {}: bottom overlap {:.3e}", start, s.name, n,
                     res.bottom_overlap_end);
      }
    }
  }
}

// a square on the ground, at its rest overlap for DEM (two corners)
auto square_at_rest(world& w, const scheme_t& s, double alpha, double vx = 0.) {
  const double y0 = is_dem(s) ? r - m * g * std::cos(alpha) / (2 * kn) : r;
  return std::pair{w.add_square(0., y0, vx), y0};
}

void incline_stick(const std::filesystem::path& dir) {
  using namespace siconos;
  auto f = open_csv(dir, "incline_stick", "scheme,h,t,x,y");
  const double alpha = 20 * deg;
  for (const auto& s : all_schemes) {
    const double h = is_dem(s) ? 1e-4 : 5e-3;
    world w(scene_for(s, 0.5, h, 1., alpha));
    auto [p, y0] = square_at_rest(w, s, alpha);
    w.finalize();
    cost c;
    const std::size_t every = std::max<std::size_t>(1, std::lround(5e-3 / h));
    while (w.simulation.has_next_event()) {
      step(w, c);
      const auto st = w.simulation.current_step();
      if (st % every == 0) {
        const auto& q = storage::attr<"q">(p, st);
        std::println(f, "{},{},{},{},{}", s.name, h, st * h, q(0), q(1) - y0);
      }
    }
    std::println(stderr, "incline_stick {}: x(T) = {:.3e}", s.name,
                 storage::attr<"q">(p, w.simulation.current_step())(0));
  }
}

// a square on the ground, run until tmax: its last x, its largest lift-off
struct square_result {
  cost c;
  bool finite = true;
  double x_end = 0;
  double max_lift = 0;
};

// path(t, x, y - y0, theta): the trajectory, if given
square_result square(const scheme_t& s, double mu, double h, double tmax, double alpha,
                     double vx,
                     std::function<void(double, double, double, double)> path = {}) {
  using namespace siconos;
  world w(scene_for(s, mu, h, tmax, alpha));
  auto [p, y0] = square_at_rest(w, s, alpha, vx);
  w.finalize();
  square_result res;
  while (w.simulation.has_next_event()) {
    step(w, res.c);
    const auto& q = storage::attr<"q">(p, w.simulation.current_step());
    if (path) path(w.simulation.current_step() * h, q(0), q(1) - y0, q(2));
    // blown up: away from the ground by more than half its size
    if (!std::isfinite(q(1)) || std::abs(q(1) - y0) > r) {
      res.finite = false;
      break;
    }
    res.max_lift = std::max(res.max_lift, q(1) - y0);
    res.x_end = q(0);
  }
  return res;
}

const std::vector<double> square_steps = {1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4};

// the trajectories of the square, one run per scheme: h = 5e-3, 1e-4 for
// DEM, every 5e-3
void square_paths(const std::filesystem::path& dir, const std::string& name, double mu,
                  double tmax, double alpha, double vx) {
  auto f = open_csv(dir, name, "scheme,h,t,x,lift,theta");
  for (const auto& s : all_schemes) {
    const double h = is_dem(s) ? 1e-4 : 5e-3;
    const std::size_t every = std::max<std::size_t>(1, std::lround(5e-3 / h));
    std::size_t k = 0;
    square(s, mu, h, tmax, alpha, vx, [&](double t, double x, double lift, double theta) {
      if (++k % every == 0) std::println(f, "{},{},{},{},{},{}", s.name, h, t, x, lift, theta);
    });
  }
}

void sliding(const std::filesystem::path& dir) {
  auto f =
      open_csv(dir, "sliding",
               "scheme,h,steps,wall_s,iterations,finite,distance,distance_exact,max_lift");
  const double v0 = 2., mu = 0.3;
  const double d_exact = v0 * v0 / (2 * mu * g);
  for (const auto& s : all_schemes) {
    for (double h : square_steps) {
      const auto res = square(s, mu, h, v0 / (mu * g) + 0.3, 0., v0);
      std::println(f, "{},{},{},{},{},{},{},{},{}", s.name, h, res.c.steps, res.c.wall(),
                   res.c.iterations, res.finite, res.x_end, d_exact, res.max_lift);
      std::println(stderr, "sliding {} h {:.0e}: distance {:.6f} (exact {:.6f}), lift {:.2e}",
                   s.name, h, res.x_end, d_exact, res.max_lift);
    }
  }
  square_paths(dir, "sliding_path", mu, v0 / (mu * g) + 0.3, 0., v0);
}

void incline_slide(const std::filesystem::path& dir) {
  auto f = open_csv(dir, "incline_slide",
                    "scheme,h,steps,wall_s,iterations,finite,x_end,x_exact,max_lift");
  const double alpha = 30 * deg, mu = 0.3, tmax = 1.;
  const double x_exact = 0.5 * g * (std::sin(alpha) - mu * std::cos(alpha)) * tmax * tmax;
  for (const auto& s : all_schemes) {
    for (double h : square_steps) {
      const auto res = square(s, mu, h, tmax, alpha, 0.);
      std::println(f, "{},{},{},{},{},{},{},{},{}", s.name, h, res.c.steps, res.c.wall(),
                   res.c.iterations, res.finite, res.x_end, x_exact, res.max_lift);
      std::println(stderr,
                   "incline_slide {} h {:.0e}: x(T) {:.6f} (exact {:.6f}), lift {:.2e}",
                   s.name, h, res.x_end, x_exact, res.max_lift);
    }
  }
  square_paths(dir, "incline_slide_path", mu, tmax, alpha, 0.);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::println(stderr, "usage: {} <output directory> [case...]", argv[0]);
    return 2;
  }
  const std::filesystem::path dir = argv[1];
  std::filesystem::create_directories(dir);
  const std::vector<std::pair<std::string, std::function<void(const std::filesystem::path&)>>>
      cases = {{"stability", stability},
               {"column", columns},
               {"incline_stick", incline_stick},
               {"sliding", sliding},
               {"incline_slide", incline_slide}};
  std::vector<std::string> wanted(argv + 2, argv + argc);
  for (const auto& [name, run] : cases) {
    if (wanted.empty() || std::find(wanted.begin(), wanted.end(), name) != wanted.end()) {
      run(dir);
    }
  }
  return 0;
}
