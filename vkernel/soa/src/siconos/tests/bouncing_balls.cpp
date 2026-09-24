
#include <chrono>
#include <numeric>
#include <fstream>
#include <print>

#include "siconos/siconos.hpp"

namespace siconos::config {
namespace some = siconos::storage::some;
struct ball : model::lagrangian_ds {};
struct lcp : simul::nonsmooth_problem<LinearComplementarityProblem> {};
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct nslaw : model::newton_impact_friction {};
struct relation : model::lagrangian_r<nslaw::size> {};
struct interaction : simul::interaction<nslaw, relation> {};
struct topo : simul::topology<ball, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct simulation : simul::time_stepping<td, osi, osnspb> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>>;
};

struct data_t
    : storage::make<
          env, simulation,
          storage::pattern::wrap<some::unbounded_collection, ball>,
          storage::pattern::wrap<some::bounded_collection, relation,
                                 some::indice_value<1>>,
          storage::pattern::wrap<some::unbounded_collection, interaction>,
          storage::with_properties<
              storage::time_invariant<storage::attr_t<ball, "fext">>,
              storage::diagonal<storage::attr_t<ball, "mass_matrix">>,
              storage::assembled_diagonal<
                  storage::attr_t<typename osi::assembled_osi_t,
                                  "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

int main(int argc, char* argv[])
{
  using namespace siconos;

  config::data_t data;

  // unsigned int nDof = 3;         // degrees of freedom for the ball
  double t0 = 0;               // initial computation time
  double tmax = 1;             // final computation time
  double h = 0.005;            // time step
  double position_init = 1.0;  // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double theta = 0.5;          // theta for MoreauJeanOSI integrator
  double radius = 0.1;         // Ball radius
  double m = 1.;               // Ball mass
  double g = 9.81;             // Gravity

  unsigned int nballs = 10;
  std::print("====> Model loading ...\n");

  // ---------------------------
  // -- The dynamical_systems --
  // ---------------------------
  for (unsigned int i = 0; i < nballs; ++i) {
    auto ball = storage::add<config::ball>(data);
    ball.q() = {position_init * (i + 1), 0, 0};
    ball.velocity() = {velocity_init, 0, 0};
    ball.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;
    ball.fext() = {-m * g, 0., 0.};
  }

  for (auto ball : storage::handles<config::ball>(data, 0)) {
    std::print("ball:{} , ball.q()={}\n", ball.index().value(), ball.q()[0]);
  }
  // ------------------
  // -- The relation --
  // ------------------

  // -- Lagrangian relation --
  auto relation_f = storage::add<config::relation>(data);
  auto relation_b = storage::add<config::relation>(data);
  relation_f.h_matrix() << 1., 0., 0., 0., 1., -radius;
  relation_b.h_matrix() << -1., 0., 0., 0., 1., -radius;
  relation_f.b()(0) = -radius;
  relation_b.b()(0) = -2 * radius;

  // -- nslaw --
  double e = 0.9;
  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;
  nslaw.mu() = 0.;

  //  auto lcp = storage::add<config::lcp>(data);
  //  lcp.create();
  auto fc2d = storage::add<config::fc2d>(data);
  fc2d.create();
  fc2d.instance()->dimension = 2;
  fc2d.instance()->mu = 0;

  // ------------------
  // --- Simulation ---
  // ------------------
  auto simulation = storage::add<config::simulation>(data);

  simulation.one_step_integrator().theta() = theta;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = t0;
  simulation.time_discretization().h() = h;

  simulation.time_discretization().tmax() = tmax;
  // -- set the formulation for the one step nonsmooth problem --
  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;

  // -- set the options --
  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_NSGS);
  osnspb.options() = so;

  auto balls = storage::handles<config::ball>(data, 0);

  auto first_ball = (balls | view::take(1)).front();
  //    views::transform([&simulation, &radius, &relation_f, &nslaw](auto
  //    first_ball)
  //    {
  auto interaction = simulation.topology().link(first_ball);
  //  interaction.h_matrix1() << 1., 0., 0., 0., 1., -radius;
  //  interaction.h_matrix2() << 1., 0., 0., 0., 1., -radius;
  interaction.relation() = relation_f;
  interaction.nslaw() = nslaw;
  //    });

  for (auto [ball1, ball2] : view::zip(balls, balls | view::drop(1))) {
    std::print("new interaction ball<->ball : {} {}\n", ball1.index().value(),
               ball2.index().value());
    auto interaction = simulation.topology().link(ball1, ball2);
    //    interaction.h_matrix1() << -1., 0., 0., 0., 1., -radius;
    //    interaction.h_matrix2() << 1., 0., 0., 0., 1., -radius;
    interaction.relation() = relation_b;
    interaction.nslaw() = nslaw;
  };

  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");
  balls = storage::handles<config::ball>(data, 0);
  auto ball1 = (balls | view::take(1)).front();
  auto ball2 = (balls | view::take(2)).back();

  // iteration matrix & h_matrices
  simulation.initialize();

  std::ofstream result_file("result-many.dat");

  std::print(result_file,
             "{:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
             simulation.current_step() * simulation.time_step(),
             storage::attr<"q">(ball1, simulation.current_step())(0),
             storage::attr<"q">(ball2, simulation.current_step())(0),
             storage::attr<"velocity">(ball1, simulation.current_step())(0),
             storage::attr<"velocity">(ball2, simulation.current_step())(0),
             0., 0.);

  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  while (simulation.has_next_event()) {
    auto ninvds = simulation.compute_one_step();

    double p01, p02, lambda1, lambda2;
    if (ninvds > 1) {
      p01 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                       0)(0);
      p02 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                       0)(1);

      lambda1 = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 0)(0);

      lambda2 = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 1)(0);
    }
    else if (ninvds == 1) {
      p01 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                       0)(0);
      p02 = 0;

      lambda1 = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 0)(0);

      lambda2 = 0;
    }
    else {
      p01 = 0;
      p02 = 0;
      lambda1 = 0;
      lambda2 = 0;
    }

    std::print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
               simulation.current_step() * simulation.time_step(),
               storage::attr<"q">(ball1, simulation.current_step())(0),
               storage::attr<"q">(ball2, simulation.current_step())(0),
               storage::attr<"velocity">(ball1, simulation.current_step())(0),
               storage::attr<"velocity">(ball2, simulation.current_step())(0),
               p01, p02, lambda1, lambda2);
  }

  std::print("Computation Time \n");
  end = std::chrono::system_clock::now();
  int elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  std::print("Computation time : {} ms \n", elapsed);

  //  io::close(fd);
}
