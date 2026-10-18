// storage::remove keeps the dynamic properties of the items consistent.
//
// Three balls, ball i at q(0) = i.  Balls 0 and 2 get a bc_velocities_0 set
// (a dynamic property attached by the topology).  Removing ball 0 moves ball 2
// into index 0 (swap and pop): its BC set must follow it, and nothing may be
// left at index 2.  Returns the number of failed checks.

#include <initializer_list>
#include <print>

#include "siconos/config/config_builder.hpp"
#include "siconos/siconos.hpp"

namespace siconos::config {
struct ball : model::lagrangian_ds {};
struct lcp : simul::nonsmooth_problem<LinearComplementarityProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<lcp> {};
struct nslaw : model::newton_impact {};
struct relation : model::lagrangian_r<nslaw::size> {};
struct interaction : simul::interaction<nslaw, relation> {};
struct topo : simul::topology<ball, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct simulation : simul::time_stepping<td, osi, osnspb> {};

struct data_t : decltype(storage()
                             .with_param<"dof", 3>()
                             .with_items<simulation>()
                             .with_time_invariant<ball, "fext">()
                             .with_diagonal<ball, "mass_matrix">()
                             .with_assembled_diagonal<osi, "mass_matrix_assembled">()
#ifdef BB_DATABASE
                             .build()
#else
                             .build_world()
#endif
                ){};
}  // namespace siconos::config

int main() {
  namespace storage = siconos::storage;
  namespace config = siconos::config;

  config::data_t data;
  int failures = 0;
  auto check = [&failures](bool ok, const char* what) {
    std::print("{} {}\n", ok ? "ok  " : "FAIL", what);
    failures += !ok;
  };

  auto ball0 = storage::add<config::ball>(data);
  auto ball1 = storage::add<config::ball>(data);
  auto ball2 = storage::add<config::ball>(data);
  ball0.q() = {0., 0., 0.};
  ball1.q() = {1., 0., 0.};
  ball2.q() = {2., 0., 0.};

  // the BC sets are Eigen vectors of indices: fill them element by element
  // (bc_t{2, 22} would be Eigen's rows x cols constructor)
  using bc_t = std::decay_t<decltype(storage::prop<"bc_velocities_0", config::topo>(ball0))>;
  auto make_bc = [](std::initializer_list<std::size_t> values) {
    bc_t bc;
    bc.resize(values.size());
    std::size_t i = 0;
    for (auto v : values) bc[i++] = v;
    return bc;
  };
  auto is_bc = [](const bc_t& bc, std::initializer_list<std::size_t> values) {
    if (static_cast<std::size_t>(bc.size()) != values.size()) return false;
    std::size_t i = 0;
    for (auto v : values) {
      if (bc[i++] != v) return false;
    }
    return true;
  };
  storage::prop<"bc_velocities_0", config::topo>(ball0) = make_bc({0});
  storage::prop<"bc_velocities_0", config::topo>(ball2) = make_bc({2, 22});

  storage::remove(data, ball0);

  auto& qs = storage::attr_values<config::ball, "q">(data, 0);
  auto& bcs = storage::prop_values<config::ball, "bc_velocities_0", config::topo>(data, 0);

  check(std::size(qs) == 2, "two balls left");
  check(qs[0](0) == 2., "ball 2 moved to index 0");
  check(qs[1](0) == 1., "ball 1 still at index 1");
  check(bcs.contains(0) && is_bc(bcs.at(0), {2, 22}),
        "ball 2's BC set followed it to index 0");
  check(!bcs.contains(1), "ball 1 still has no BC set");
  check(!bcs.contains(2), "nothing left at index 2");

  std::print("{} failure(s)\n", failures);
  return failures;
}
