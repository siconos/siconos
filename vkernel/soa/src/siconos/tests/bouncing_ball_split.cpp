#include "bouncing_ball_split.hpp"

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <fstream>
#include <print>

#include "ref_check.hpp"

namespace {

class BouncingBallSplitTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(BouncingBallSplitTest);
  CPPUNIT_TEST(testBouncingBallSplit);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testBouncingBallSplit()
  {
    namespace storage = siconos::storage;
    namespace config = siconos::config;
    using siconos::storage::handle;
    using siconos::storage::pattern::param_val;

    config::data_t data;

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

    std::print("====> Model loading ...\n");

    // --------------------------
    // -- The dynamical_system --
    // --------------------------
    handle ball = storage::add<config::ball>(data);

    ball.q() = {position_init, 0, 0};
    ball.velocity() = {velocity_init, 0, 0};
    ball.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

    // -- Set external forces (weight) --
    ball.fext() = {-m * g, 0., 0.};

    // ------------------
    // -- The relation --
    // ------------------

    // -- Lagrangian relation --
    handle relation = storage::add<config::relation>(data);
    relation.h_matrix() = {-1.0, 0., 0.};

    // -- nslaw --
    double e = 0.9;
    handle nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = e;

    handle lcp = storage::add<config::lcp>(data);
    lcp.create(SICONOS_LCP_LEMKE);

    // ------------------
    // --- Simulation ---
    // ------------------
    handle simulation = storage::add<config::simulation>(data);

    simulation.one_step_integrator().theta() = theta;
    simulation.one_step_integrator().constraint_activation_threshold() = 0.;
    simulation.time_discretization().t0() = t0;
    simulation.time_discretization().h() = h;

    simulation.time_discretization().tmax() = tmax;
    // -- set the formulation for the one step nonsmooth problem --
    handle osnspb = simulation.one_step_nonsmooth_problem();
    osnspb.problem() = lcp;

    // -- set the options --
    handle so = storage::add<siconos::simul::solver_options>(data);
    so.create();
    osnspb.options() = so;

    // Interaction ball-floor
    handle interaction = simulation.topology().link(ball);
    // interaction.h_matrix1() = {1., 0., 0.};

    simulation.topology().bind(interaction, relation);
    interaction.nslaw() = nslaw;

    // =========================== End of model definition
    // ===========================
    // ================================= Computation
    // =================================

    //  auto fd = io::open<ascii>("bouncing_ball_split.dat");

    // fix this for constant fext
    simulation_initialize(data, simulation.index().value());

    std::ofstream result_file("bouncing_ball_split.dat");

    std::print(result_file, "{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
               simulation.current_step() * simulation.time_step(),
               storage::attr<"q">(ball, simulation.current_step())(0),
               storage::attr<"velocity">(ball, simulation.current_step())(0),
               0., 0.);

    while (simulation_has_next_event(data, simulation.index().value())) {
      unsigned int ninvds =
          simulation_compute_one_step(data, simulation.index().value());

      double p0, lambda;
      if (ninvds > 0) {
        p0 = get_vector(
            simulation.one_step_integrator().p0_vector_assembled(), 0)(0);
        lambda = get_vector(
            simulation.one_step_integrator().lambda_vector_assembled(), 0)(0);
      }
      else {
        p0 = 0;
        lambda = 0;
      }

      std::print(
          result_file, "{:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
          simulation.current_step() * simulation.time_step(),
          storage::attr<"q">(ball, simulation.current_step())(0),
          storage::attr<"velocity">(ball, simulation.current_step())(0), p0,
          lambda);
    }
    //  io::close(fd);

    result_file.close();
    CPPUNIT_ASSERT_REF(siconos::test::read_ref("bouncing_ball_split.dat"),
                       "bouncing_ball_split.ref", 1e-8);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(BouncingBallSplitTest);

}  // namespace
