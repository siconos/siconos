// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "TimeSteppingSchemesTest.hpp"

#include <cmath>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include "BoundaryCondition.hpp"
#include "FirstOrderType2R.hpp"
#include "LagrangianSparseLinearTIDS.hpp"
#include "SiconosKernel.hpp"
#include "io.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION(TimeSteppingSchemesTest);

namespace {

using siconos::algebra::SiconosDenseMatrix;
using siconos::algebra::SiconosVector;
using siconos::modeling::Interaction;
using siconos::modeling::NonSmoothDynamicalSystem;
using siconos::simulation::TimeDiscretisation;

namespace mod = siconos::modeling;
namespace integ = siconos::integrators;
namespace sim = siconos::simulation;
namespace osns = siconos::nonsmooth_formulations;

constexpr double g = 9.81;
constexpr int number_of_samples = 10;
// Relative tolerance of the comparison, column by column:
// max_i |x_ij - ref_ij| <= tolerance * (1 + max_i |ref_ij|).
constexpr double tolerance = 1e-8;

// ---------------------------------------------------------------------------------------
// A test case: the system, the simulation, the storage the relations alias, and what to
// record at each step.
// ---------------------------------------------------------------------------------------
struct Case {
  std::shared_ptr<NonSmoothDynamicalSystem> nsds;
  std::shared_ptr<sim::Simulation> simulation;
  std::vector<std::shared_ptr<mod::DynamicalSystem>> dss;
  std::vector<std::shared_ptr<Interaction>> interactions;
  // LagrangianLinearTIR and friends alias their matrices/vectors: keep them alive here.
  std::deque<siconos::algebra::SiconosMatrix> matrices;
  std::deque<SiconosVector> vectors;
  // level of lambda to record (1 for the velocity-level schemes, 0 for Schatzman-Paoli)
  unsigned int lambda_level{1};
  // additional values to record (e.g. for interactions created during the simulation)
  std::function<void(std::vector<double>&)> extra;
  // record the multipliers (not when they are not unique, e.g. a box resting on 4 points)
  bool record_lambda{true};
};

void append(std::vector<double>& row, const Eigen::Ref<const Eigen::VectorXd>& x) {
  for (Eigen::Index i = 0; i < x.size(); ++i) row.push_back(x(i));
}

/* One row of the trajectory: time, then q and v of each dynamical system (x for a
   first-order one), then for each interaction the first component of y[0] (the gap) and
   of lambda[level] (the normal multiplier). */
std::vector<double> record(Case& c) {
  std::vector<double> row{c.simulation->nextTime()};
  for (auto& ds : c.dss) {
    // Lagrangian (dense or sparse) and Newton-Euler systems (velocity = twist)
    if (auto sods = std::dynamic_pointer_cast<mod::SecondOrderDS>(ds)) {
      append(row, sods->q_read());
      append(row, sods->velocity_read());
    } else if (auto fods = std::dynamic_pointer_cast<mod::FirstOrderNonLinearDS>(ds)) {
      append(row, fods->x_read());
    }
  }
  for (auto& inter : c.interactions) {
    row.push_back((*inter->y(0))(0));
    if (!c.record_lambda) continue;
    auto lambda = inter->lambda(c.lambda_level);
    row.push_back((lambda && lambda->size() > 0) ? (*lambda)(0) : 0.0);
  }
  if (c.extra) c.extra(row);
  return row;
}

/* Run nsteps steps and return 10 rows of the trajectory, evenly spaced, the last one at
   the final time. With SICONOS_TEST_DUMP_TRAJECTORIES set, the whole trajectory is also
   written to <name>.dat in the working directory (for inspection). */
SiconosDenseMatrix run(Case& c, int nsteps, const std::string& name = "") {
  std::vector<std::vector<double>> rows;
  try {
    for (int k = 0; k < nsteps; ++k) {
      c.simulation->advanceToEvent();
      rows.push_back(record(c));
      c.simulation->processEvents();
    }
  } catch (...) {
    std::cerr << "\n" << name << ": exception at step " << rows.size() << "\n";
    siconos::exception::process();  // prints the details of a siconos exception
    throw;
  }
  if (std::getenv("SICONOS_TEST_DUMP_TRAJECTORIES") && !name.empty()) {
    SiconosDenseMatrix all(static_cast<Eigen::Index>(rows.size()),
                           static_cast<Eigen::Index>(rows[0].size()));
    for (std::size_t k = 0; k < rows.size(); ++k)
      for (std::size_t j = 0; j < rows[k].size(); ++j)
        all(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(j)) = rows[k][j];
    siconos::algebra::io::write(name + ".dat", all, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
  }
  SiconosDenseMatrix samples(number_of_samples, static_cast<Eigen::Index>(rows[0].size()));
  for (int s = 0; s < number_of_samples; ++s) {
    auto k = static_cast<std::size_t>(((s + 1) * nsteps) / number_of_samples - 1);
    for (std::size_t j = 0; j < rows[k].size(); ++j)
      samples(s, static_cast<Eigen::Index>(j)) = rows[k][j];
  }
  return samples;
}

/* Compare the samples with <name>.ref (copied in the working directory of the test), or
   (re)generate the reference in the source directory if SICONOS_REGENERATE_TEST_REFS is
   set. */
void check(const std::string& name, const SiconosDenseMatrix& samples) {
  const std::string file = name + ".ref";
  CPPUNIT_ASSERT_MESSAGE(name + ": NaN in the trajectory", samples.allFinite());
  if (std::getenv("SICONOS_REGENERATE_TEST_REFS")) {
    std::string src = std::string(SICONOS_TEST_REF_DIR) + "/" + file;
    siconos::algebra::io::write(src, samples);
    siconos::algebra::io::write(file, samples);
    std::cout << "\n  [regenerated " << src << "]";
    return;
  }
  if (!std::filesystem::exists(file))
    CPPUNIT_FAIL(name + ": no reference file " + file +
                 " (run the test with SICONOS_REGENERATE_TEST_REFS=1 to create it)");
  auto ref = siconos::algebra::io::readDenseMatrix(file);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": number of samples", ref.rows(), samples.rows());
  CPPUNIT_ASSERT_EQUAL_MESSAGE(name + ": number of columns", ref.cols(), samples.cols());
  for (Eigen::Index j = 0; j < ref.cols(); ++j) {
    double scale = 1.0 + ref.col(j).cwiseAbs().maxCoeff();
    Eigen::Index i;
    double err = (samples.col(j) - ref.col(j)).cwiseAbs().maxCoeff(&i) / scale;
    if (!(err <= tolerance)) {
      std::ostringstream msg;
      msg.precision(16);
      msg << name << ": column " << j << ", sample " << i << " (t = " << samples(i, 0)
          << "): " << samples(i, j) << " instead of " << ref(i, j) << " (relative error "
          << err << " > " << tolerance << ")";
      CPPUNIT_FAIL(msg.str());
    }
  }
}

void setTolerance(osns::OneStepNSProblem& pb, double tol, int maxiter = 0) {
  pb.numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = tol;
  if (maxiter > 0) pb.numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = maxiter;
}

// ---------------------------------------------------------------------------------------
// Simulation setups
// ---------------------------------------------------------------------------------------
enum class Scheme {
  moreau_jean,
  combined_projection,
  combined_projection_from_origin,
  direct_projection,
  d1minus_acceleration,
  d1minus_velocity,
  schatzman_paoli,
};

/* Time stepping with one LCP (or MLCP if `mixed`) at the velocity level, and the
   position-level problem of the projection schemes. */
/* `tune`, if given, configures the velocity-level problem (storage, assembly). */
void setupSimulation(Case& c, Scheme scheme, double h, double theta = 0.5, bool mixed = false,
                     std::function<void(osns::LinearOSNS&)> tune = nullptr) {
  auto td = std::make_shared<TimeDiscretisation>(c.nsds->t0(), h);
  auto velocity_problem = [&]() -> std::shared_ptr<osns::OneStepNSProblem> {
    std::shared_ptr<osns::LinearOSNS> pb;
    if (mixed) {
      pb = std::make_shared<osns::MLCP>(SICONOS_MLCP_ENUM);
      setTolerance(*pb, 1e-12);
    } else {
      pb = std::make_shared<osns::LCP>();
    }
    if (tune) tune(*pb);
    return pb;
  };
  auto position_problem = [&]() {
    auto pb = std::make_shared<osns::MLCPProjectOnConstraints>(SICONOS_MLCP_ENUM);
    // the enumerative solver accepts z, w >= -tol: keep tol far below the gaps
    setTolerance(*pb, 1e-12, 2000);
    return pb;
  };
  switch (scheme) {
    case Scheme::moreau_jean: {
      auto osi = std::make_shared<integ::MoreauJeanOSI>(theta);
      c.simulation = std::make_shared<sim::TimeStepping>(c.nsds, td, osi, velocity_problem());
      break;
    }
    case Scheme::combined_projection:
    case Scheme::combined_projection_from_origin: {
      auto osi = std::make_shared<integ::MoreauJeanCombinedProjectionOSI>(theta);
      auto s = std::make_shared<sim::TimeSteppingCombinedProjection>(
          c.nsds, td, osi, velocity_problem(), position_problem());
      s->setConstraintTol(1e-10);
      s->setConstraintTolUnilateral(1e-10);
      if (scheme == Scheme::combined_projection_from_origin)
        s->setProjectionUpdate(sim::ProjectionUpdate::from_origin);
      c.simulation = s;
      break;
    }
    case Scheme::direct_projection: {
      auto osi = std::make_shared<integ::MoreauJeanDirectProjectionOSI>(theta);
      auto s = std::make_shared<sim::TimeSteppingDirectProjection>(
          c.nsds, td, osi, velocity_problem(), position_problem(), 0);
      s->setConstraintTol(1e-10);
      s->setConstraintTolUnilateral(1e-10);
      c.simulation = s;
      break;
    }
    case Scheme::d1minus_acceleration:
    case Scheme::d1minus_velocity: {
      auto type = scheme == Scheme::d1minus_acceleration
                      ? integ::D1MinusLinearOSI::Type::halfexplicit_acceleration_level
                      : integ::D1MinusLinearOSI::Type::halfexplicit_velocity_level;
      auto osi = std::make_shared<integ::D1MinusLinearOSI>(type);
      auto s = std::make_shared<sim::TimeSteppingD1Minus>(c.nsds, td, 2);
      s->insertIntegrator(osi);
      s->insertNonSmoothProblem(std::make_shared<osns::LCP>(), sim::SICONOS_OSNSP_TS_VELOCITY);
      s->insertNonSmoothProblem(std::make_shared<osns::LCP>(),
                                sim::SICONOS_OSNSP_TS_VELOCITY + 1);
      c.simulation = s;
      break;
    }
    case Scheme::schatzman_paoli: {
      auto osi = std::make_shared<integ::SchatzmanPaoliOSI>(theta);
      c.simulation =
          std::make_shared<sim::TimeStepping>(c.nsds, td, osi, std::make_shared<osns::LCP>());
      c.lambda_level = 0;
      break;
    }
  }
  if (auto ts = std::dynamic_pointer_cast<sim::TimeStepping>(c.simulation)) {
    ts->setNewtonTolerance(1e-12);
    ts->setNewtonMaxIteration(50);
  }
}

// ---------------------------------------------------------------------------------------
// Systems
// ---------------------------------------------------------------------------------------

constexpr double beads_radius = 0.1;
constexpr double beads_e = 0.8;

/* Interaction between two beads i < i+1: y = z_{i+1} - z_i - 2R. */
std::shared_ptr<Interaction> beadsInteraction(Case& c) {
  auto& H = c.matrices.emplace_back(1, 2);
  H << -1.0, 1.0;
  auto& b = c.vectors.emplace_back(1);
  b << -2 * beads_radius;
  return std::make_shared<Interaction>(std::make_shared<mod::NewtonImpactNSL>(beads_e),
                                       std::make_shared<mod::LagrangianLinearTIR>(H, b));
}

/* Interaction manager creating the bead/bead interactions when the beads get closer than
   `link_gap`, and removing them when they get farther than `unlink_gap`, as a contact
   detection would. */
class BeadsContactDetection : public sim::InteractionManager {
  Case& case_;
  double link_gap_, unlink_gap_;

 public:
  std::vector<std::shared_ptr<Interaction>> pairs;  // pairs[i]: beads i and i+1 (or null)
  int links{0}, unlinks{0};

  BeadsContactDetection(Case& c, double link_gap, double unlink_gap)
      : case_{c}, link_gap_{link_gap}, unlink_gap_{unlink_gap}, pairs(c.dss.size() - 1) {}

  void updateInteractions(std::shared_ptr<sim::Simulation> simulation) override {
    for (std::size_t i = 0; i + 1 < case_.dss.size(); ++i) {
      double zi = std::static_pointer_cast<mod::LagrangianDS>(case_.dss[i])->q_read()(0);
      double zj = std::static_pointer_cast<mod::LagrangianDS>(case_.dss[i + 1])->q_read()(0);
      double gap = zj - zi - 2 * beads_radius;
      if (!pairs[i] && gap < link_gap_) {
        pairs[i] = beadsInteraction(case_);
        simulation->link(pairs[i], case_.dss[i], case_.dss[i + 1]);
        links++;
      } else if (pairs[i] && gap > unlink_gap_) {
        simulation->unlink(pairs[i]);
        pairs[i].reset();
        unlinks++;
      }
    }
  }
};

/* Column of 3 beads (1 dof each, radius 0.1) above the ground, falling with different
   velocities: impacts with the ground and between beads (simultaneous ones included),
   then resting contacts. 3 unilateral contacts, linear relations, a single restitution
   coefficient (as in the hots-python models, see validate_with_hots.py).
   If `managed`, only the ground contact is given: the bead/bead interactions are created
   and removed during the simulation by BeadsContactDetection (see beadsManaged). */
Case beads(bool managed = false) {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  const double R = beads_radius;
  const double z0[3] = {0.15, 0.45, 0.8};
  const double v0[3] = {0.0, -1.0, -2.0};
  for (int i = 0; i < 3; ++i) {
    SiconosVector q{1}, v{1}, f{1};
    q << z0[i];
    v << v0[i];
    f << -g;
    siconos::algebra::SiconosDenseMatrix M =
        siconos::algebra::SiconosDenseMatrix::Identity(1, 1);
    auto ds = std::make_shared<mod::LagrangianLinearTIDS>(q, v, M, siconos::algebra::copy_t);
    ds->setConstantFext(f, siconos::algebra::copy_t);
    c.nsds->insertDynamicalSystem(ds);
    c.dss.push_back(ds);
  }
  // ground: y = z0 - R
  auto& H0 = c.matrices.emplace_back(1, 1);
  H0 << 1.0;
  auto& b0 = c.vectors.emplace_back(1);
  b0 << -R;
  auto inter0 =
      std::make_shared<Interaction>(std::make_shared<mod::NewtonImpactNSL>(beads_e),
                                    std::make_shared<mod::LagrangianLinearTIR>(H0, b0));
  c.nsds->link(inter0, c.dss[0]);
  c.interactions.push_back(inter0);
  if (managed) return c;
  // bead/bead: y = z_{i+1} - z_i - 2R
  for (int i = 0; i < 2; ++i) {
    auto inter = beadsInteraction(c);
    c.nsds->link(inter, c.dss[i], c.dss[i + 1]);
    c.interactions.push_back(inter);
  }
  return c;
}

/* Chain of 5 masses on a line, linked by springs and dampers (sparse stiffness and damping
   matrices, LagrangianSparseLinearTIDS), moving towards two obstacles: the last mass hits
   a wall, the middle one a stop behind it. The first mass is either clamped by a boundary
   condition (`clamped`) or attached to a fixed wall by a spring. */
Case springChain(bool clamped, bool dense = false) {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  const int n = 5;
  const double k = 100.0, damping = 0.01;
  SiconosVector q = SiconosVector::Zero(n), v = SiconosVector::Constant(n, 1.0);
  if (clamped) v(0) = 0.0;
  std::vector<Eigen::Triplet<double>> tk, tm;
  for (int i = 0; i < n; ++i) {
    tm.emplace_back(i, i, 1.0);
    // spring between i-1 (or the wall for i = 0) and i
    tk.emplace_back(i, i, k);
    if (i > 0) {
      tk.emplace_back(i - 1, i - 1, k);
      tk.emplace_back(i - 1, i, -k);
      tk.emplace_back(i, i - 1, -k);
    }
  }
  siconos::algebra::SiconosSparseMatrix M(n, n), K(n, n);
  M.setFromTriplets(tm.begin(), tm.end());
  K.setFromTriplets(tk.begin(), tk.end());
  siconos::algebra::SiconosSparseMatrix C = damping * K;
  std::shared_ptr<mod::SecondOrderDS> ds;
  if (dense) {  // for the integrators restricted to LagrangianLinearTIDS (SchatzmanPaoliOSI)
    auto dds = std::make_shared<mod::LagrangianLinearTIDS>(
        q, v, siconos::algebra::SiconosDenseMatrix(M), siconos::algebra::copy_t);
    dds->setStiffnessMatrix(siconos::algebra::SiconosDenseMatrix(K), siconos::algebra::copy_t);
    dds->setDampingMatrix(siconos::algebra::SiconosDenseMatrix(C), siconos::algebra::copy_t);
    ds = dds;
  } else {
    auto sds =
        std::make_shared<mod::LagrangianSparseLinearTIDS>(q, v, M, siconos::algebra::copy_t);
    sds->setStiffnessMatrix(K, siconos::algebra::copy_t);
    sds->setDampingMatrix(C, siconos::algebra::copy_t);
    ds = sds;
  }
  if (clamped)
    ds->setBoundaryConditions(
        std::make_shared<mod::BoundaryCondition>(mod::BoundaryCondition::Indices{0}));
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  // wall at q_4 = 0.05 (y = 0.05 - q_4), stop at q_2 = 0.08 (y = 0.08 - q_2)
  for (auto [dof, position] : {std::pair{4, 0.05}, std::pair{2, 0.08}}) {
    auto& H = c.matrices.emplace_back(1, n);
    H.setZero();
    H(0, dof) = -1.0;
    auto& b = c.vectors.emplace_back(1);
    b << position;
    auto inter =
        std::make_shared<Interaction>(std::make_shared<mod::NewtonImpactNSL>(0.5),
                                      std::make_shared<mod::LagrangianLinearTIR>(H, b));
    c.nsds->link(inter, ds);
    c.interactions.push_back(inter);
  }
  return c;
}

/* Two independent linear oscillators (diagonal stiffness) hitting obstacles. */
Case oscillators() {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  SiconosVector q{2}, v{2}, k{2};
  q << 0.5, 0.3;
  v << 0.0, 1.0;
  k << 100.0, 400.0;
  auto ds =
      std::make_shared<mod::LagrangianLinearDiagonalDS>(q, v, k, siconos::algebra::copy_t);
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  for (int i = 0; i < 2; ++i) {
    // obstacle at q_i = -0.2: y = q_i + 0.2
    auto& H = c.matrices.emplace_back(1, 2);
    H.setZero();
    H(0, i) = 1.0;
    auto& b = c.vectors.emplace_back(1);
    b << 0.2;
    auto inter =
        std::make_shared<Interaction>(std::make_shared<mod::NewtonImpactNSL>(0.7),
                                      std::make_shared<mod::LagrangianLinearTIR>(H, b));
    c.nsds->link(inter, ds);
    c.interactions.push_back(inter);
  }
  return c;
}

/* Disk (z, x, theta) thrown with a horizontal velocity onto the ground, Coulomb friction
   in 2D: bounces, sliding, then rolling without sliding. */
Case disk2D() {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  const double R = 0.1, m = 1.0;
  SiconosVector q{3}, v{3}, f{3};
  q << R + 0.2, 0.0, 0.0;
  v << 0.0, 2.0, 0.0;
  f << -m * g, 0.0, 0.0;
  SiconosDenseMatrix M = SiconosDenseMatrix::Zero(3, 3);
  M(0, 0) = m;
  M(1, 1) = m;
  M(2, 2) = 0.5 * m * R * R;
  auto ds = std::make_shared<mod::LagrangianLinearTIDS>(q, v, M, siconos::algebra::copy_t);
  ds->setConstantFext(f, siconos::algebra::copy_t);
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  // normal: z - R, tangent: velocity of the contact point x' + R theta'
  auto& H = c.matrices.emplace_back(2, 3);
  H << 1.0, 0.0, 0.0, 0.0, 1.0, R;
  auto& b = c.vectors.emplace_back(2);
  b << -R, 0.0;
  auto inter = std::make_shared<Interaction>(
      std::make_shared<mod::NewtonImpactFrictionNSL>(0.5, 0.0, 0.3, 2),
      std::make_shared<mod::LagrangianLinearTIR>(H, b));
  c.nsds->link(inter, ds);
  c.interactions.push_back(inter);
  return c;
}

/* Sphere (x, y, z, rotation vector) with initial sliding velocity and spin, Coulomb
   friction in 3D. */
Case sphere3D() {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  const double R = 0.1, m = 1.0, I = 0.4 * m * R * R;
  SiconosVector q{6}, v{6}, f{6};
  q << 0.0, 0.0, R + 0.3, 0.0, 0.0, 0.0;
  v << 1.5, 0.5, 0.0, 10.0, 0.0, 3.0;
  f << 0.0, 0.0, -m * g, 0.0, 0.0, 0.0;
  SiconosDenseMatrix M = SiconosDenseMatrix::Zero(6, 6);
  M.diagonal() << m, m, m, I, I, I;
  auto ds = std::make_shared<mod::LagrangianLinearTIDS>(q, v, M, siconos::algebra::copy_t);
  ds->setConstantFext(f, siconos::algebra::copy_t);
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  // contact point at (0, 0, -R) from the center: u_c = v + w x r = v + (-R w_y, R w_x, 0)
  auto& H = c.matrices.emplace_back(3, 6);
  H << 0, 0, 1, 0, 0, 0,  //
      1, 0, 0, 0, -R, 0,  //
      0, 1, 0, R, 0, 0;
  auto& b = c.vectors.emplace_back(3);
  b << -R, 0.0, 0.0;
  auto inter = std::make_shared<Interaction>(
      std::make_shared<mod::NewtonImpactFrictionNSL>(0.6, 0.0, 0.3, 3),
      std::make_shared<mod::LagrangianLinearTIR>(H, b));
  c.nsds->link(inter, ds);
  c.interactions.push_back(inter);
  return c;
}

/* Rocking block (x, y, theta) released with an angle on one corner: rotation about the
   corner, impacts alternating between the two corners. Nonlinear relations. */
Case rockingBlock() {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  constexpr double L = 0.2, Hb = 0.1, m = 1.0;
  const double theta0 = 0.3;
  SiconosVector q{3}, v{3}, f{3};
  // corner 1 (-L/2, -H/2 in the block frame) on the ground
  q << 0.5 * L * std::cos(theta0) - 0.5 * Hb * std::sin(theta0),
      0.5 * L * std::sin(theta0) + 0.5 * Hb * std::cos(theta0), theta0;
  v << 0.0, 0.0, 0.0;
  f << 0.0, -m * g, 0.0;
  SiconosDenseMatrix M = SiconosDenseMatrix::Zero(3, 3);
  M.diagonal() << m, m, m * (L * L + Hb * Hb) / 12.0;
  auto ds = std::make_shared<mod::LagrangianLinearTIDS>(q, v, M, siconos::algebra::copy_t);
  ds->setConstantFext(f, siconos::algebra::copy_t);
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  for (double s : {-1.0, 1.0}) {
    // corner (s L/2, -H/2): y = q1 + s L/2 sin(theta) - H/2 cos(theta)
    auto r = std::make_shared<mod::LagrangianScleronomousR>();
    r->setComputehFunction(
        [s](const siconos::algebra::BlockVector& q, Eigen::Ref<SiconosVector> y) {
          y(0) = q(1) + s * 0.5 * L * std::sin(q(2)) - 0.5 * Hb * std::cos(q(2));
        });
    r->setComputeContactMatrixFunction(
        [s](const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::MapType> H) {
          H(0, 0) = 0.0;
          H(0, 1) = 1.0;
          H(0, 2) = s * 0.5 * L * std::cos(q(2)) + 0.5 * Hb * std::sin(q(2));
        });
    r->setComputeContactMatrixDotFunction([s](const siconos::algebra::BlockVector& q,
                                              const siconos::algebra::BlockVector& qdot,
                                              Eigen::Ref<siconos::algebra::MapType> Hdot) {
      Hdot(0, 0) = 0.0;
      Hdot(0, 1) = 0.0;
      Hdot(0, 2) = (-s * 0.5 * L * std::sin(q(2)) + 0.5 * Hb * std::cos(q(2))) * qdot(2);
    });
    auto inter = std::make_shared<Interaction>(std::make_shared<mod::NewtonImpactNSL>(0.5), r);
    c.nsds->link(inter, ds);
    c.interactions.push_back(inter);
  }
  return c;
}

/* Double pendulum (angles from the vertical, nonlinear mass, gyroscopic and gravity
   forces) released against a wall x = 0 hit by both masses (e = 0.8 and e = 0). */
Case doublePendulum() {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  constexpr double m1 = 1.0, m2 = 1.0, l1 = 1.0, l2 = 1.0;
  SiconosVector q{2}, v{2};
  q << 0.6, 1.2;
  v << 0.0, 0.0;
  auto ds = std::make_shared<mod::LagrangianDS>(q, v, siconos::algebra::copy_t);
  ds->setComputeMassFunction(
      [](const Eigen::Ref<const SiconosVector>& q, Eigen::Ref<siconos::algebra::MapType> M) {
        M(0, 0) = (m1 + m2) * l1 * l1;
        M(0, 1) = M(1, 0) = m2 * l1 * l2 * std::cos(q(0) - q(1));
        M(1, 1) = m2 * l2 * l2;
      });
  // M(q) q'' + fgyr(q, v) + fint(q) = p,
  // fgyr = m2 l1 l2 sin(q0 - q1) (v1^2, -v0^2), fint = g (( m1 + m2) l1 sin q0, m2 l2 sin q1)
  ds->setComputeFgyrFunction([](const Eigen::Ref<const SiconosVector>& v,
                                const Eigen::Ref<const SiconosVector>& q,
                                Eigen::Ref<siconos::algebra::MapVectorType> fgyr) {
    double c = m2 * l1 * l2 * std::sin(q(0) - q(1));
    fgyr(0) = c * v(1) * v(1);
    fgyr(1) = -c * v(0) * v(0);
  });
  ds->setComputeJacobianFgyrOver_qFunction([](const Eigen::Ref<const SiconosVector>& v,
                                              const Eigen::Ref<const SiconosVector>& q,
                                              Eigen::Ref<siconos::algebra::MapType> J) {
    double c = m2 * l1 * l2 * std::cos(q(0) - q(1));
    J(0, 0) = c * v(1) * v(1);
    J(0, 1) = -c * v(1) * v(1);
    J(1, 0) = -c * v(0) * v(0);
    J(1, 1) = c * v(0) * v(0);
  });
  ds->setComputeJacobianFgyrOver_velocityFunction([](const Eigen::Ref<const SiconosVector>& v,
                                                     const Eigen::Ref<const SiconosVector>& q,
                                                     Eigen::Ref<siconos::algebra::MapType> J) {
    double c = m2 * l1 * l2 * std::sin(q(0) - q(1));
    J(0, 0) = 0.0;
    J(0, 1) = 2 * c * v(1);
    J(1, 0) = -2 * c * v(0);
    J(1, 1) = 0.0;
  });
  ds->setComputeFintFunction([](const Eigen::Ref<const SiconosVector>&,
                                const Eigen::Ref<const SiconosVector>& q, double,
                                Eigen::Ref<siconos::algebra::MapVectorType> fint) {
    fint(0) = (m1 + m2) * g * l1 * std::sin(q(0));
    fint(1) = m2 * g * l2 * std::sin(q(1));
  });
  ds->setComputeJacobianFintOver_qFunction([](const Eigen::Ref<const SiconosVector>&,
                                              const Eigen::Ref<const SiconosVector>& q, double,
                                              Eigen::Ref<siconos::algebra::MapType> J) {
    J(0, 0) = (m1 + m2) * g * l1 * std::cos(q(0));
    J(0, 1) = 0.0;
    J(1, 0) = 0.0;
    J(1, 1) = m2 * g * l2 * std::cos(q(1));
  });
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  // wall x = 0 (x = l1 sin q0 for the first mass, + l2 sin q1 for the second)
  for (int k = 0; k < 2; ++k) {
    auto r = std::make_shared<mod::LagrangianScleronomousR>();
    r->setComputehFunction(
        [k](const siconos::algebra::BlockVector& q, Eigen::Ref<SiconosVector> y) {
          y(0) = l1 * std::sin(q(0)) + (k == 1 ? l2 * std::sin(q(1)) : 0.0);
        });
    r->setComputeContactMatrixFunction(
        [k](const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::MapType> H) {
          H(0, 0) = l1 * std::cos(q(0));
          H(0, 1) = k == 1 ? l2 * std::cos(q(1)) : 0.0;
        });
    r->setComputeContactMatrixDotFunction([k](const siconos::algebra::BlockVector& q,
                                              const siconos::algebra::BlockVector& qdot,
                                              Eigen::Ref<siconos::algebra::MapType> Hdot) {
      Hdot(0, 0) = -l1 * std::sin(q(0)) * qdot(0);
      Hdot(0, 1) = k == 1 ? -l2 * std::sin(q(1)) * qdot(1) : 0.0;
    });
    auto inter = std::make_shared<Interaction>(
        std::make_shared<mod::NewtonImpactNSL>(k == 0 ? 0.8 : 0.0), r);
    c.nsds->link(inter, ds);
    c.interactions.push_back(inter);
  }
  return c;
}

/* Point mass in Cartesian coordinates, attached by a rigid rod of length L (bilateral
   constraint x^2 + y^2 - L^2 = 0), swinging down onto a floor y = -0.8 L (unilateral). */
Case bilateralPendulum() {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  constexpr double L = 1.0, m = 1.0;
  SiconosVector q{2}, v{2}, f{2};
  q << L, 0.0;
  v << 0.0, 0.0;
  f << 0.0, -m * g;
  SiconosDenseMatrix M = m * SiconosDenseMatrix::Identity(2, 2);
  auto ds = std::make_shared<mod::LagrangianLinearTIDS>(q, v, M, siconos::algebra::copy_t);
  ds->setConstantFext(f, siconos::algebra::copy_t);
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  auto rod = std::make_shared<mod::LagrangianScleronomousR>();
  rod->setComputehFunction(
      [](const siconos::algebra::BlockVector& q, Eigen::Ref<SiconosVector> y) {
        y(0) = q(0) * q(0) + q(1) * q(1) - L * L;
      });
  rod->setComputeContactMatrixFunction(
      [](const siconos::algebra::BlockVector& q, Eigen::Ref<siconos::algebra::MapType> H) {
        H(0, 0) = 2 * q(0);
        H(0, 1) = 2 * q(1);
      });
  auto inter_rod =
      std::make_shared<Interaction>(std::make_shared<mod::EqualityConditionNSL>(1), rod);
  c.nsds->link(inter_rod, ds);
  c.interactions.push_back(inter_rod);
  auto& H = c.matrices.emplace_back(1, 2);
  H << 0.0, 1.0;
  auto& b = c.vectors.emplace_back(1);
  b << 0.8 * L;
  auto inter_floor =
      std::make_shared<Interaction>(std::make_shared<mod::NewtonImpactNSL>(0.5),
                                    std::make_shared<mod::LagrangianLinearTIR>(H, b));
  c.nsds->link(inter_floor, ds);
  c.interactions.push_back(inter_floor);
  return c;
}

/* Contact of a point of a Newton-Euler body (given in the body frame) with the plane
   z = 0. Base is NewtonEulerNormalVelocityR (frictionless) or NewtonEulerVelocityR (with
   Coulomb friction). */
template <class Base>
class PointOnPlaneR : public Base {
  siconos::algebra::SiconosVector3 point_;

 public:
  explicit PointOnPlaneR(const siconos::algebra::SiconosVector3& point)
      : Base{}, point_{point} {}

  void computeh(const Eigen::Ref<const siconos::algebra::SiconosVector7>& q,
                const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>&,
                Eigen::Ref<SiconosVector> y) override {
    Eigen::Quaterniond rotation(q(3), q(4), q(5), q(6));
    Eigen::Vector3d p = q.head<3>() + rotation.normalized() * point_;
    this->contactPoint1_ = p;
    this->contactPoint2_ << p(0), p(1), 0.0;
    this->nc_ << 0.0, 0.0, 1.0;
    y.setZero();
    y(0) = p(2);
  }
};

/* Cube (side 0.2) dropped with a tilt, a horizontal velocity and a spin onto the plane
   z = 0: impacts of one or several corners, rocking, sliding/sticking, rest. One
   contact per corner (8 interactions). */
Case cube(bool friction) {
  Case c;
  c.record_lambda = false;  // 4 contacts at rest on a rigid face: not unique
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  const double a = 0.2, m = 1.0;
  siconos::algebra::SiconosVector7 q;
  Eigen::Quaterniond tilt(Eigen::AngleAxisd(0.3, Eigen::Vector3d(1, 1, 0).normalized()));
  q << 0.0, 0.0, 0.25, tilt.w(), tilt.x(), tilt.y(), tilt.z();
  siconos::algebra::SiconosVector6 v;
  v << 0.5, 0.0, 0.0, 0.0, 2.0, 1.0;
  siconos::algebra::SiconosMatrix33 I =
      (m * a * a / 6.0) * siconos::algebra::SiconosMatrix33::Identity();
  auto ds = std::make_shared<mod::NewtonEulerDS>(q, v, m, I, siconos::algebra::copy_t);
  siconos::algebra::SiconosVector3 weight(0.0, 0.0, -m * g);
  ds->setConstantFext(weight, siconos::algebra::copy_t);
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  for (double sx : {-1.0, 1.0})
    for (double sy : {-1.0, 1.0})
      for (double sz : {-1.0, 1.0}) {
        siconos::algebra::SiconosVector3 corner(sx * a / 2, sy * a / 2, sz * a / 2);
        std::shared_ptr<mod::Relation> r;
        std::shared_ptr<mod::NonSmoothLaw> nslaw;
        if (friction) {
          r = std::make_shared<PointOnPlaneR<mod::NewtonEulerVelocityR>>(corner);
          nslaw = std::make_shared<mod::NewtonImpactFrictionNSL>(0.5, 0.0, 0.3, 3);
        } else {
          r = std::make_shared<PointOnPlaneR<mod::NewtonEulerNormalVelocityR>>(corner);
          nslaw = std::make_shared<mod::NewtonImpactNSL>(0.5);
        }
        auto inter = std::make_shared<Interaction>(nslaw, r);
        c.nsds->link(inter, ds);
        c.interactions.push_back(inter);
      }
  return c;
}

/* Mass (x, y, z) held in the plane by springs and dampers (stiffness and damping
   matrices), launched on a plane with Coulomb friction: stick-slip oscillations. Dense
   (LagrangianLinearTIDS) or sparse (LagrangianSparseLinearTIDS) storage. */
Case frictionOscillator(bool sparse) {
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  const double m = 1.0, kx = 50.0, ky = 80.0, cd = 0.2;
  SiconosVector q{3}, v{3}, f{3};
  q << 0.1, -0.05, 0.02;
  v << 0.0, 1.0, 0.0;
  f << 0.0, 0.0, -m * g;
  SiconosDenseMatrix M = m * SiconosDenseMatrix::Identity(3, 3),
                     K = SiconosDenseMatrix::Zero(3, 3), C = SiconosDenseMatrix::Zero(3, 3);
  K(0, 0) = kx;
  K(1, 1) = ky;
  K(0, 1) = K(1, 0) = 10.0;
  C(0, 0) = C(1, 1) = cd;
  std::shared_ptr<mod::LagrangianDS> lds;
  if (sparse) {
    siconos::algebra::SiconosSparseMatrix Ms = M.sparseView(), Ks = K.sparseView(),
                                          Cs = C.sparseView();
    auto ds =
        std::make_shared<mod::LagrangianSparseLinearTIDS>(q, v, Ms, siconos::algebra::copy_t);
    ds->setStiffnessMatrix(Ks, siconos::algebra::copy_t);
    ds->setDampingMatrix(Cs, siconos::algebra::copy_t);
    ds->setConstantFext(f, siconos::algebra::copy_t);
    c.dss.push_back(ds);
  } else {
    auto ds = std::make_shared<mod::LagrangianLinearTIDS>(q, v, M, siconos::algebra::copy_t);
    ds->setStiffnessMatrix(K, siconos::algebra::copy_t);
    ds->setDampingMatrix(C, siconos::algebra::copy_t);
    ds->setConstantFext(f, siconos::algebra::copy_t);
    c.dss.push_back(ds);
  }
  c.nsds->insertDynamicalSystem(c.dss[0]);
  auto& H = c.matrices.emplace_back(3, 3);
  H << 0, 0, 1,  //
      1, 0, 0,   //
      0, 1, 0;
  auto inter = std::make_shared<Interaction>(
      std::make_shared<mod::NewtonImpactFrictionNSL>(0.0, 0.0, 0.3, 3),
      std::make_shared<mod::LagrangianLinearTIR>(H));
  c.nsds->link(inter, c.dss[0]);
  c.interactions.push_back(inter);
  return c;
}

/* Contact of a point of body 1 (given in its frame) with a face of body 2 (point and
   outward normal of the face given in the frame of body 2): the face is taken as an
   infinite plane. */
template <class Base>
class PointOnFaceR : public Base {
  siconos::algebra::SiconosVector3 point_, face_point_, face_normal_;

 public:
  PointOnFaceR(const siconos::algebra::SiconosVector3& point,
               const siconos::algebra::SiconosVector3& face_point,
               const siconos::algebra::SiconosVector3& face_normal)
      : Base{}, point_{point}, face_point_{face_point}, face_normal_{face_normal} {}

  void computeh(const Eigen::Ref<const siconos::algebra::SiconosVector7>& q1,
                const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector7>>& q2,
                Eigen::Ref<SiconosVector> y) override {
    assert(q2);
    Eigen::Quaterniond r1(q1(3), q1(4), q1(5), q1(6)),
        r2((*q2)(3), (*q2)(4), (*q2)(5), (*q2)(6));
    Eigen::Vector3d p = q1.head<3>() + r1.normalized() * point_;
    Eigen::Vector3d a = q2->head<3>() + r2.normalized() * face_point_;
    Eigen::Vector3d n = r2.normalized() * face_normal_;
    double gap = (p - a).dot(n);
    this->contactPoint1_ = p;
    this->contactPoint2_ = p - gap * n;
    this->nc_ = n;
    y.setZero();
    y(0) = gap;
  }
};

/* Two boxes (side 0.2): the lower one lies on the ground (3 bottom points), the upper one
   is dropped with a tilt and a spin on the top face of the lower one (3 bottom points):
   3 interactions between the same pair of bodies, the lower box is pushed into the
   ground.
   \warning With TimeSteppingCombinedProjection, the simulation diverges when the upper box
   settles on its 3 contacts (t = 0.164: the upper box is projected to z = 1.3e7), see
   setupSimulation(c, Scheme::combined_projection, 4e-3) on this system. (Not used with
   friction: resting on 4 frictional contacts is an indeterminate problem on which the friction
   solvers converge very slowly.) */
Case stack(bool friction) {
  Case c;
  c.record_lambda = false;  // 4 contacts at rest on a rigid face: not unique
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 1.0);
  const double a = 0.2, m = 1.0;
  siconos::algebra::SiconosMatrix33 I =
      (m * a * a / 6.0) * siconos::algebra::SiconosMatrix33::Identity();
  siconos::algebra::SiconosVector3 weight(0.0, 0.0, -m * g);
  Eigen::Quaterniond tilt(Eigen::AngleAxisd(0.15, Eigen::Vector3d(1, 2, 0).normalized()));
  siconos::algebra::SiconosVector7 q_low, q_up;
  q_low << 0.0, 0.0, a / 2, 1.0, 0.0, 0.0, 0.0;
  q_up << 0.01, -0.01, 1.5 * a + 0.05, tilt.w(), tilt.x(), tilt.y(), tilt.z();
  siconos::algebra::SiconosVector6 v_low = siconos::algebra::SiconosVector6::Zero(), v_up;
  v_up << 0.1, 0.0, -0.5, 0.0, 0.0, 1.5;
  auto low =
      std::make_shared<mod::NewtonEulerDS>(q_low, v_low, m, I, siconos::algebra::copy_t);
  auto up = std::make_shared<mod::NewtonEulerDS>(q_up, v_up, m, I, siconos::algebra::copy_t);
  for (auto& ds : {low, up}) {
    ds->setConstantFext(weight, siconos::algebra::copy_t);
    c.nsds->insertDynamicalSystem(ds);
    c.dss.push_back(ds);
  }
  auto nslaw = [&]() -> std::shared_ptr<mod::NonSmoothLaw> {
    if (friction) return std::make_shared<mod::NewtonImpactFrictionNSL>(0.3, 0.0, 0.4, 3);
    return std::make_shared<mod::NewtonImpactNSL>(0.3);
  };
  // 3 points of the bottom face of each box (2 corners and the middle of the opposite
  // edge, around the center of mass): with the 4 corners, the 4 coplanar contacts are
  // redundant and the projection of the combined projection scheme (G G^T singular)
  // blows up when the upper box lands.
  for (auto [sx, sy] : {std::pair{-1.0, -1.0}, std::pair{1.0, -1.0}, std::pair{0.0, 1.0}}) {
    siconos::algebra::SiconosVector3 corner(sx * a / 2, sy * a / 2, -a / 2);
    std::shared_ptr<mod::Relation> ground, face;
    siconos::algebra::SiconosVector3 top(0.0, 0.0, a / 2), ez(0.0, 0.0, 1.0);
    if (friction) {
      ground = std::make_shared<PointOnPlaneR<mod::NewtonEulerVelocityR>>(corner);
      face = std::make_shared<PointOnFaceR<mod::NewtonEulerVelocityR>>(corner, top, ez);
    } else {
      ground = std::make_shared<PointOnPlaneR<mod::NewtonEulerNormalVelocityR>>(corner);
      face = std::make_shared<PointOnFaceR<mod::NewtonEulerNormalVelocityR>>(corner, top, ez);
    }
    auto inter_ground = std::make_shared<Interaction>(nslaw(), ground);
    c.nsds->link(inter_ground, low);
    c.interactions.push_back(inter_ground);
    auto inter_face = std::make_shared<Interaction>(nslaw(), face);
    c.nsds->link(inter_face, up, low);
    c.interactions.push_back(inter_face);
  }
  return c;
}

/* Time stepping with a (global) frictional contact problem, solved to convergence. */
void setupFrictionSimulation(Case& c, int dim, bool global, double h, double theta = 0.5,
                             std::function<void(osns::LinearOSNS&)> tune = nullptr) {
  auto td = std::make_shared<TimeDiscretisation>(c.nsds->t0(), h);
  std::shared_ptr<osns::OneStepNSProblem> pb;
  std::shared_ptr<integ::OneStepIntegrator> osi;
  if (global) {
    pb =
        std::make_shared<osns::GlobalFrictionContact>(dim, SICONOS_GLOBAL_FRICTION_3D_NSGS_WR);
    osi = std::make_shared<integ::MoreauJeanGOSI>(theta);
  } else {
    pb = std::make_shared<osns::FrictionContact>(
        dim, dim == 2 ? SICONOS_FRICTION_2D_LEMKE : SICONOS_FRICTION_3D_NSGS);
    osi = std::make_shared<integ::MoreauJeanOSI>(theta);
  }
  // the results must not depend on where an iterative solver stops
  if (tune) tune(static_cast<osns::LinearOSNS&>(*pb));
  setTolerance(*pb, 1e-14, 100000);
  auto ts = std::make_shared<sim::TimeStepping>(c.nsds, td, osi, pb);
  ts->setNewtonTolerance(1e-12);
  ts->setNewtonMaxIteration(50);
  c.simulation = ts;
}

}  // namespace

// ---------------------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------------------

void TimeSteppingSchemesTest::beadsMoreauJean() {
  auto c = beads();
  setupSimulation(c, Scheme::moreau_jean, 2e-3);
  check("beadsMoreauJean", run(c, 300, "beadsMoreauJean"));
}

void TimeSteppingSchemesTest::beadsCombinedProjection() {
  auto c = beads();
  setupSimulation(c, Scheme::combined_projection, 2e-3);
  check("beadsCombinedProjection", run(c, 300, "beadsCombinedProjection"));
}

void TimeSteppingSchemesTest::beadsDirectProjection() {
  auto c = beads();
  setupSimulation(c, Scheme::direct_projection, 2e-3);
  check("beadsDirectProjection", run(c, 300, "beadsDirectProjection"));
}

void TimeSteppingSchemesTest::beadsD1MinusAccelerationLevel() {
  auto c = beads();
  setupSimulation(c, Scheme::d1minus_acceleration, 2e-3);
  check("beadsD1MinusAccelerationLevel", run(c, 300, "beadsD1MinusAccelerationLevel"));
}

void TimeSteppingSchemesTest::beadsD1MinusVelocityLevel() {
  auto c = beads();
  setupSimulation(c, Scheme::d1minus_velocity, 2e-3);
  check("beadsD1MinusVelocityLevel", run(c, 300, "beadsD1MinusVelocityLevel"));
}

void TimeSteppingSchemesTest::beadsSchatzmanPaoli() {
  auto c = beads();
  setupSimulation(c, Scheme::schatzman_paoli, 2e-3);
  check("beadsSchatzmanPaoli", run(c, 300, "beadsSchatzmanPaoli"));
}

void TimeSteppingSchemesTest::oscillatorsMoreauJeanBilbao() {
  auto c = oscillators();
  auto td = std::make_shared<TimeDiscretisation>(0.0, 2e-3);
  c.simulation = std::make_shared<sim::TimeStepping>(
      c.nsds, td, std::make_shared<integ::MoreauJeanBilbaoOSI>(),
      std::make_shared<osns::LCP>());
  check("oscillatorsMoreauJeanBilbao", run(c, 400, "oscillatorsMoreauJeanBilbao"));
}

void TimeSteppingSchemesTest::disk2DFrictionMoreauJean() {
  auto c = disk2D();
  setupFrictionSimulation(c, 2, false, 5e-3);
  check("disk2DFrictionMoreauJean", run(c, 200, "disk2DFrictionMoreauJean"));
}

void TimeSteppingSchemesTest::sphere3DFrictionMoreauJean() {
  auto c = sphere3D();
  setupFrictionSimulation(c, 3, false, 5e-3);
  check("sphere3DFrictionMoreauJean", run(c, 200, "sphere3DFrictionMoreauJean"));
}

void TimeSteppingSchemesTest::sphere3DFrictionMoreauJeanGOSI() {
  auto c = sphere3D();
  setupFrictionSimulation(c, 3, true, 5e-3);
  check("sphere3DFrictionMoreauJeanGOSI", run(c, 200, "sphere3DFrictionMoreauJeanGOSI"));
}

void TimeSteppingSchemesTest::rockingBlockMoreauJean() {
  auto c = rockingBlock();
  setupSimulation(c, Scheme::moreau_jean, 5e-3);
  check("rockingBlockMoreauJean", run(c, 200, "rockingBlockMoreauJean"));
}

void TimeSteppingSchemesTest::rockingBlockCombinedProjection() {
  auto c = rockingBlock();
  setupSimulation(c, Scheme::combined_projection, 5e-3);
  check("rockingBlockCombinedProjection", run(c, 200, "rockingBlockCombinedProjection"));
}

void TimeSteppingSchemesTest::rockingBlockCombinedProjectionFromOrigin() {
  auto c = rockingBlock();
  setupSimulation(c, Scheme::combined_projection_from_origin, 5e-3);
  check("rockingBlockCombinedProjectionFromOrigin",
        run(c, 200, "rockingBlockCombinedProjectionFromOrigin"));
}

void TimeSteppingSchemesTest::rockingBlockDirectProjection() {
  auto c = rockingBlock();
  setupSimulation(c, Scheme::direct_projection, 5e-3);
  check("rockingBlockDirectProjection", run(c, 200, "rockingBlockDirectProjection"));
}

void TimeSteppingSchemesTest::rockingBlockD1Minus() {
  auto c = rockingBlock();
  setupSimulation(c, Scheme::d1minus_acceleration, 5e-3);
  check("rockingBlockD1Minus", run(c, 200, "rockingBlockD1Minus"));
}

void TimeSteppingSchemesTest::doublePendulumMoreauJean() {
  auto c = doublePendulum();
  setupSimulation(c, Scheme::moreau_jean, 2e-3);
  check("doublePendulumMoreauJean", run(c, 400, "doublePendulumMoreauJean"));
}

void TimeSteppingSchemesTest::doublePendulumD1Minus() {
  auto c = doublePendulum();
  setupSimulation(c, Scheme::d1minus_acceleration, 2e-3);
  check("doublePendulumD1Minus", run(c, 400, "doublePendulumD1Minus"));
}

void TimeSteppingSchemesTest::bilateralPendulumMoreauJean() {
  auto c = bilateralPendulum();
  setupSimulation(c, Scheme::moreau_jean, 2e-3, 0.5, true);
  check("bilateralPendulumMoreauJean", run(c, 400, "bilateralPendulumMoreauJean"));
}

void TimeSteppingSchemesTest::bilateralPendulumCombinedProjection() {
  auto c = bilateralPendulum();
  setupSimulation(c, Scheme::combined_projection, 2e-3, 0.5, true);
  check("bilateralPendulumCombinedProjection",
        run(c, 400, "bilateralPendulumCombinedProjection"));
}

void TimeSteppingSchemesTest::cubeFrictionMoreauJean() {
  // until t = 0.5: afterwards the cube rests on 4 frictional contacts, an indeterminate
  // problem on which NSGS converges very slowly (and the Newton solvers are much slower)
  auto c = cube(true);
  setupFrictionSimulation(c, 3, false, 4e-3);
  check("cubeFrictionMoreauJean", run(c, 125, "cubeFrictionMoreauJean"));
}

void TimeSteppingSchemesTest::cubeFrictionMoreauJeanGOSI() {
  // until t = 0.5: afterwards the cube rests on 4 frictional contacts, an indeterminate
  // problem on which NSGS converges very slowly (and the Newton solvers are much slower)
  auto c = cube(true);
  setupFrictionSimulation(c, 3, true, 4e-3);
  check("cubeFrictionMoreauJeanGOSI", run(c, 125, "cubeFrictionMoreauJeanGOSI"));
}

void TimeSteppingSchemesTest::cubeMoreauJean() {
  auto c = cube(false);
  setupSimulation(c, Scheme::moreau_jean, 4e-3);
  check("cubeMoreauJean", run(c, 200, "cubeMoreauJean"));
}

void TimeSteppingSchemesTest::cubeDirectProjection() {
  auto c = cube(false);
  setupSimulation(c, Scheme::direct_projection, 4e-3);
  check("cubeDirectProjection", run(c, 200, "cubeDirectProjection"));
}

void TimeSteppingSchemesTest::cubeCombinedProjection() {
  auto c = cube(false);
  setupSimulation(c, Scheme::combined_projection, 4e-3);
  check("cubeCombinedProjection", run(c, 200, "cubeCombinedProjection"));
}

void TimeSteppingSchemesTest::freeRigidBodyMoreauJean() {
  // rotation close to the intermediate principal axis: unstable (Dzhanibekov effect)
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 2.0);
  siconos::algebra::SiconosVector7 q;
  q << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
  siconos::algebra::SiconosVector6 v;
  // growth rate of the instability w sqrt((I2 - I1)(I3 - I2) / (I1 I3)) = 2.9 / s: the
  // perturbation 0.1 leads to a flip around t = 1.35
  v << 0.0, 0.0, 0.0, 0.1, 5.0, 0.1;
  siconos::algebra::SiconosMatrix33 I = siconos::algebra::SiconosMatrix33::Zero();
  I.diagonal() << 1.0, 2.0, 3.0;
  auto ds = std::make_shared<mod::NewtonEulerDS>(q, v, 1.0, I, siconos::algebra::copy_t);
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  setupSimulation(c, Scheme::moreau_jean, 1e-2);
  check("freeRigidBodyMoreauJean", run(c, 200, "freeRigidBodyMoreauJean"));
}

void TimeSteppingSchemesTest::circuitRLCDEulerMoreau() { circuitRLCD(false); }

void TimeSteppingSchemesTest::circuitRLCDEulerMoreauResiduals() {
  // Newton loop with the residual-based convergence criteria on y and r
  circuitRLCD(true);
}

void TimeSteppingSchemesTest::circuitRLCD(bool residuals) {
  // RLC circuit with an ideal diode (complementarity), as in CircuitRLCD of the tutorials
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 5e-3);
  const double L = 1e-2, C = 1e-6, R = 1e3;
  SiconosVector x0{2};
  x0 << 10.0, 0.0;
  auto ds = std::make_shared<mod::FirstOrderLinearDS>(x0, siconos::algebra::copy_t);
  auto& A = c.matrices.emplace_back(2, 2);
  A << 0.0, -1.0 / C, 1.0 / L, 0.0;
  ds->setConstantA(A, siconos::algebra::alias_t);
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);
  auto& Cm = c.matrices.emplace_back(1, 2);
  Cm << -1.0, 0.0;
  auto& Bm = c.matrices.emplace_back(2, 1);
  Bm << -1.0 / C, 0.0;
  auto& Dm = c.matrices.emplace_back(1, 1);
  Dm << R;
  auto r = std::make_shared<mod::FirstOrderLinearTIR>(Cm, Bm);
  r->setConstantD(Dm);
  auto inter =
      std::make_shared<Interaction>(std::make_shared<mod::ComplementarityConditionNSL>(1), r);
  c.nsds->link(inter, ds);
  c.interactions.push_back(inter);
  auto td = std::make_shared<TimeDiscretisation>(0.0, 1e-5);
  c.simulation = std::make_shared<sim::TimeStepping>(
      c.nsds, td, std::make_shared<integ::EulerMoreauOSI>(0.5), std::make_shared<osns::LCP>());
  c.lambda_level = 0;
  std::string name = residuals ? "circuitRLCDEulerMoreauResiduals" : "circuitRLCDEulerMoreau";
  if (residuals) {
    auto ts = std::static_pointer_cast<sim::TimeStepping>(c.simulation);
    ts->setComputeResiduY(true);
    ts->setComputeResiduR(true);
  }
  check(name, run(c, 500, name));
}

void TimeSteppingSchemesTest::stepSystemEulerMoreauType2R() {
  // Piecewise-linear gene regulatory network (SignSystem of the tutorials, biology/
  // StepSystem): x' = f(x) + g(lambda), y = h(x), relay law -1 <= lambda <= 1 (sign of y),
  // with a nonlinear g (FirstOrderType2R relation). EulerMoreauOSI with the Newton loop
  // on the linearised relation and the residual of the output. The trajectory reaches the
  // switching surfaces and slides on them.
  Case c;
  c.nsds = std::make_shared<NonSmoothDynamicalSystem>(0.0, 3.0);
  SiconosVector x0{2};
  x0 << 1.0, 6.0;
  auto ds = std::make_shared<mod::FirstOrderNonLinearDS>(x0, siconos::algebra::copy_t);
  ds->setComputefVectorFunction([](const Eigen::Ref<const SiconosVector>& x, double,
                                   Eigen::Ref<siconos::algebra::MapVectorType> f) {
    f(0) = -4.5 * x(0);
    f(1) = -1.5 * x(1);
  });
  ds->setComputeJacobianfOver_xFunction([](const Eigen::Ref<const SiconosVector>&, double,
                                           Eigen::Ref<siconos::algebra::MapType> J) {
    J.setZero();
    J(0, 0) = -4.5;
    J(1, 1) = -1.5;
  });
  c.nsds->insertDynamicalSystem(ds);
  c.dss.push_back(ds);

  auto r = std::make_shared<mod::FirstOrderType2R>();
  r->setComputehFunction([](const siconos::algebra::BlockVector& x,
                            const Eigen::Ref<const SiconosVector>&,
                            Eigen::Ref<SiconosVector> y) {
    y(0) = x(0) - 4;
    y(1) = x(1) - 4;
    y(2) = x(0) - 8;
    y(3) = x(1) - 8;
  });
  r->setComputegFunction(
      [](const Eigen::Ref<const SiconosVector>& l, siconos::algebra::BlockVector& g) {
        g(0) = 10.0 * (1 + l(2)) * (1 - l(1));
        g(1) = 10.0 * (1 - l(0)) * (1 + l(3));
      });
  auto& C = c.matrices.emplace_back(4, 2);
  C << 1, 0, 0, 1, 1, 0, 0, 1;
  r->setConstantJacobianhOver_state(C);
  r->setComputeJacobiangOver_lambdaFunction(
      [](const Eigen::Ref<const SiconosVector>& l, Eigen::Ref<siconos::algebra::MapType> B) {
        B.setZero();
        B(1, 0) = -10.0 * (1 + l(3));
        B(0, 1) = -10.0 * (1 + l(2));
        B(0, 2) = 10.0 * (1 - l(1));
        B(1, 3) = 10.0 * (1 - l(0));
      });
  auto inter = std::make_shared<Interaction>(std::make_shared<mod::RelayNSL>(4, -1.0, 1.0), r);
  c.nsds->link(inter, ds);
  c.interactions.push_back(inter);

  auto td = std::make_shared<TimeDiscretisation>(0.0, 1e-3);
  auto osnsp = std::make_shared<osns::Relay>(SICONOS_RELAY_LEMKE);
  osnsp->numericsSolverOptions()->dparam[0] = 1e-8;
  auto ts = std::make_shared<sim::TimeStepping>(
      c.nsds, td, std::make_shared<integ::EulerMoreauOSI>(0.5), osnsp);
  ts->setComputeResiduY(true);
  ts->setUseRelativeConvergenceCriteron(false);
  ts->setNewtonTolerance(5e-4);
  ts->setNewtonMaxIteration(20);
  ts->setResetAllLambda(false);
  c.simulation = ts;
  c.lambda_level = 0;
  // the number of Newton iterations is part of the result: with a correct linearisation
  // the Newton loop converges in a few iterations (it reached the maximum, 20, at each
  // step when the linearisation was wrong)
  unsigned int maxNewton = 0;
  c.extra = [&maxNewton, ts, inter](std::vector<double>& row) {
    append(row, *inter->lambda(0));  // the 4 multipliers
    maxNewton = std::max(maxNewton, ts->getNewtonNbIterations());
    row.push_back(ts->getNewtonNbIterations());
  };
  auto samples = run(c, 3000, "stepSystemEulerMoreauType2R");
  CPPUNIT_ASSERT_MESSAGE("Newton loop not converged (max iterations reached)", maxNewton < 20);
  check("stepSystemEulerMoreauType2R", samples);
}

// ---------------------------------------------------------------------------------------
// Variants added to cover more of the time-stepping code (see the coverage report)
// ---------------------------------------------------------------------------------------

void TimeSteppingSchemesTest::beadsManagedMoreauJean() {
  // Same beads as beadsMoreauJean, but the bead/bead interactions are created and removed
  // by an interaction manager: same trajectory (the interactions exist whenever they can
  // become active).
  auto c = beads(true);
  setupSimulation(c, Scheme::moreau_jean, 2e-3);
  auto detection = std::make_shared<BeadsContactDetection>(c, 0.05, 0.1);
  c.simulation->insertInteractionManager(detection);
  c.extra = [detection](std::vector<double>& row) {
    for (auto& inter : detection->pairs) {
      row.push_back(inter ? (*inter->y(0))(0) : 0.0);
      row.push_back(inter ? (*inter->lambda(1))(0) : 0.0);
    }
  };
  auto samples = run(c, 300, "beadsManagedMoreauJean");
  CPPUNIT_ASSERT_MESSAGE("no interaction created", detection->links > 0);
  CPPUNIT_ASSERT_MESSAGE("no interaction removed", detection->unlinks > 0);
  check("beadsManagedMoreauJean", samples);
}

void TimeSteppingSchemesTest::beadsMoreauJeanSparseBlockStorage() {
  auto c = beads();
  setupSimulation(c, Scheme::moreau_jean, 2e-3, 0.5, false, [](osns::LinearOSNS& pb) {
    // the sparse block storage needs a block solver (not Lemke)
    pb.setMStorageType(NM_SPARSE_BLOCK);
    pb.setSolverId(SICONOS_LCP_NSGS_SBM);
    setTolerance(pb, 1e-14, 100000);
  });
  check("beadsMoreauJeanSparseBlockStorage", run(c, 300, "beadsMoreauJeanSparseBlockStorage"));
}

void TimeSteppingSchemesTest::sphere3DFrictionMoreauJeanReducedDirect() {
  // the reduced direct assembly is implemented for the sparse storage only (and the LCP
  // driver handles dense matrices only: the frictional contact problem is used here)
  auto c = sphere3D();
  setupFrictionSimulation(c, 3, false, 5e-3, 0.5, [](osns::LinearOSNS& pb) {
    pb.setMStorageType(NM_SPARSE);
    pb.setAssemblyType(osns::LinearOSNSAssemblyType::REDUCED_DIRECT);
  });
  check("sphere3DFrictionMoreauJeanReducedDirect",
        run(c, 200, "sphere3DFrictionMoreauJeanReducedDirect"));
}

void TimeSteppingSchemesTest::springChainClampedMoreauJean() {
  auto c = springChain(true);
  setupSimulation(c, Scheme::moreau_jean, 2e-3);
  check("springChainClampedMoreauJean", run(c, 400, "springChainClampedMoreauJean"));
}

void TimeSteppingSchemesTest::springChainSchatzmanPaoli() {
  // stiffness and damping terms of the Schatzman-Paoli scheme
  auto c = springChain(false, true);
  setupSimulation(c, Scheme::schatzman_paoli, 2.5e-3);
  check("springChainSchatzmanPaoli", run(c, 400, "springChainSchatzmanPaoli"));
}

void TimeSteppingSchemesTest::springChainMoreauJean() {
  auto c = springChain(false);
  setupSimulation(c, Scheme::moreau_jean, 2e-3);
  check("springChainMoreauJean", run(c, 400, "springChainMoreauJean"));
}

void TimeSteppingSchemesTest::frictionOscillatorMoreauJean() {
  auto c = frictionOscillator(false);
  setupFrictionSimulation(c, 3, false, 5e-3);
  check("frictionOscillatorMoreauJean", run(c, 200, "frictionOscillatorMoreauJean"));
}

void TimeSteppingSchemesTest::frictionOscillatorMoreauJeanGOSI() {
  auto c = frictionOscillator(false);
  setupFrictionSimulation(c, 3, true, 5e-3);
  check("frictionOscillatorMoreauJeanGOSI", run(c, 200, "frictionOscillatorMoreauJeanGOSI"));
}

void TimeSteppingSchemesTest::frictionOscillatorSparseMoreauJean() {
  auto c = frictionOscillator(true);
  setupFrictionSimulation(c, 3, false, 5e-3);
  check("frictionOscillatorSparseMoreauJean",
        run(c, 200, "frictionOscillatorSparseMoreauJean"));
}

void TimeSteppingSchemesTest::frictionOscillatorSparseMoreauJeanGOSI() {
  auto c = frictionOscillator(true);
  setupFrictionSimulation(c, 3, true, 5e-3);
  check("frictionOscillatorSparseMoreauJeanGOSI",
        run(c, 200, "frictionOscillatorSparseMoreauJeanGOSI"));
}

void TimeSteppingSchemesTest::sphere3DFrictionMoreauJeanGOSIGlobalReduced() {
  auto c = sphere3D();
  setupFrictionSimulation(c, 3, true, 5e-3, 0.5, [](osns::LinearOSNS& pb) {
    pb.setAssemblyType(osns::LinearOSNSAssemblyType::GLOBAL_REDUCED);
  });
  check("sphere3DFrictionMoreauJeanGOSIGlobalReduced",
        run(c, 200, "sphere3DFrictionMoreauJeanGOSIGlobalReduced"));
}

void TimeSteppingSchemesTest::bilateralPendulumDirectProjection() {
  auto c = bilateralPendulum();
  setupSimulation(c, Scheme::direct_projection, 2e-3, 0.5, true);
  check("bilateralPendulumDirectProjection", run(c, 400, "bilateralPendulumDirectProjection"));
}

void TimeSteppingSchemesTest::doublePendulumMoreauJeanNonlinearFull() {
  // Newton loop with the update of the interactions and index sets at each iteration
  // (the residual-based criteria are not implemented for MoreauJeanOSI, see
  // circuitRLCDEulerMoreauResiduals)
  auto c = doublePendulum();
  setupSimulation(c, Scheme::moreau_jean, 2e-3);
  auto ts = std::static_pointer_cast<sim::TimeStepping>(c.simulation);
  ts->setNewtonOptions(sim::TimeSteppingType::NONLINEAR_FULL);
  check("doublePendulumMoreauJeanNonlinearFull",
        run(c, 400, "doublePendulumMoreauJeanNonlinearFull"));
}

void TimeSteppingSchemesTest::stackMoreauJean() {
  auto c = stack(false);
  setupSimulation(c, Scheme::moreau_jean, 4e-3);
  check("stackMoreauJean", run(c, 200, "stackMoreauJean"));
}

void TimeSteppingSchemesTest::stackMoreauJeanSymmetricGraph() {
  // same as stackMoreauJean, with the symmetric version of the interaction graphs
  // (OneStepNSProblem::updateInteractionBlocks, symmetric case): same trajectory
  auto c = stack(false);
  c.nsds->setSymmetric(true);
  setupSimulation(c, Scheme::moreau_jean, 4e-3);
  check("stackMoreauJeanSymmetricGraph", run(c, 200, "stackMoreauJeanSymmetricGraph"));
}

void TimeSteppingSchemesTest::stackDirectProjection() {
  auto c = stack(false);
  setupSimulation(c, Scheme::direct_projection, 4e-3);
  check("stackDirectProjection", run(c, 200, "stackDirectProjection"));
}
