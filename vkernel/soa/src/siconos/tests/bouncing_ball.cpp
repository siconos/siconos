
#include "siconos/siconos.hpp"
#include "siconos/utils/print.hpp"

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

using params = map<iparam<"dof", 3>>;

struct make
    : storage::make<
          standard_environment<params>, simulation,
          storage::with_properties<
              storage::time_invariant<storage::attr_t<config::ball, "fext">>,
              storage::diagonal<storage::attr_t<config::ball, "mass_matrix">>,
              storage::assembled_diagonal<
                  storage::attr_t<typename config::osi::assembled_osi_t,
                                  "mass_matrix_assembled">>>> {};

}  // namespace siconos::config

int main(int argc, char* argv[])
{
  using namespace siconos;
  auto data = config::make();

  // unsigned int nDof = 3;         // degrees of freedom for the ball
  double t0 = 0;               // initial computation time
  double tmax = 10;            // final computation time
  double h = 0.005;            // time step
  double position_init = 1.0;  // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double theta = 0.5;          // theta for MoreauJeanOSI integrator
  double radius = 0.1;         // Ball radius
  double m = 1.;               // Ball mass
  double g = 9.81;             // Gravity

  print("====> Model loading ...\n");

  // --------------------------
  // -- The dynamical_system --
  // --------------------------
  auto ball = storage::add<config::ball>(data);

  ball.q() = {position_init, 0, 0};
  ball.velocity() = {velocity_init, 0, 0};
  ball.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  // -- Set external forces (weight) --
  ball.fext() = {-m * g, 0., 0.};

  // ------------------
  // -- The relation --
  // ------------------

  // -- Lagrangian relation --
  auto relation = storage::add<config::relation>(data);
  relation.h_matrix() = {-1.0, 0., 0.};

  // -- nslaw --
  double e = 0.9;
  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;

  auto lcp = storage::add<config::lcp>(data);
  lcp.create(SICONOS_LCP_LEMKE);

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
  osnspb.problem() = lcp;

  // -- set the options --
  auto so = storage::add<simul::solver_options>(data);
  so.create();
  osnspb.options() = so;

  // Interaction ball-floor
  auto interaction = simulation.topology().link(ball);
  // interaction.h_matrix1() = {1., 0., 0.};

  interaction.relation() = relation;
  interaction.nslaw() = nslaw;

  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");

  // fix this for constant fext
  simulation.initialize();

  auto out = fmt::output_file("result.dat");

  out.print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
            simulation.current_step() * simulation.time_step(),
            storage::attr<"q">(ball, simulation.current_step())(0),
            storage::attr<"velocity">(ball, simulation.current_step())(0), 0.,
            0.);

  while (simulation.has_next_event()) {
    auto ninvds = simulation.compute_one_step();

    double p0, lambda;
    if (ninvds > 0) {
      p0 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                      0)(0);
      lambda = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 0)(0);
    }
    else {
      p0 = 0;
      lambda = 0;
    }

    out.print("{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
              simulation.current_step() * simulation.time_step(),
              storage::attr<"q">(ball, simulation.current_step())(0),
              storage::attr<"velocity">(ball, simulation.current_step())(0),
              p0, lambda);
  }
  //  io::close(fd);
}
