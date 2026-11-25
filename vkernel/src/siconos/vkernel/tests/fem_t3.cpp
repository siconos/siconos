// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

// A linear elastic square (T3 elements, square_200.msh) falls onto a fixed
// ground: rt_rt interactions on its 8 bottom nodes, Moreau-Jean, theta = 1,
// e = 0, h = 1e-5, 1000 steps.  The same model as the AOS kernel's
// siconos-tutorials examples/mechanics/FEM/native/T3.cpp, whose output is
// T3_square_200.ref.
//
// The two trajectories agree to rounding over the whole run, the impact
// included (relative error 1.3e-14).  They used to decorrelate after the
// impact: that was not chaos but the Dirichlet conditions of the vkernel,
// whose blocked velocities were zeroed after the unconstrained solve instead
// of a solve with the operator of the free dofs (algebra::constrain_dofs).
// fem_t3.ref locks the trajectory against regression.

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <FENode.hpp>
#include <FemTools.hpp>
#include <FiniteElementModel.hpp>
#include <Material.hpp>
#include <Mesh.hpp>
#include <MeshUtils.hpp>
#include <SiconosKernel.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <print>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "ref_check.hpp"
#include "siconos/vkernel/config/config.hpp"
#include "siconos/vkernel/config/environment.hpp"
#include "siconos/vkernel/model/fem.hpp"
#include "siconos/vkernel/model/lagrangian_ds.hpp"
#include "siconos/vkernel/model/lagrangian_r.hpp"
#include "siconos/vkernel/model/nslaws.hpp"
#include "siconos/vkernel/simul/interaction.hpp"
#include "siconos/vkernel/simul/one_step_integrator.hpp"
#include "siconos/vkernel/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/vkernel/simul/time_discretization.hpp"
#include "siconos/vkernel/simul/time_stepping.hpp"
#include "siconos/vkernel/simul/topology.hpp"
#include "siconos/vkernel/storage/handle.hpp"

namespace siconos::vkernel::config {
struct fem_ds : model::elastic_lagrangian_ds {};
using ball = model::lagrangian_ds;
struct lcp : simul::nonsmooth_problem<LinearComplementarityProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<lcp> {};
using nslaw = model::newton_impact;
struct relation : model::lagrangian_r<nslaw::size> {};
struct rt_relation : model::rt_lagrangian_r {};
struct interaction : simul::interaction<nslaw, relation> {};
struct rt_rt_interaction : simul::rt_rt_interaction<nslaw, rt_relation> {};
struct topo : simul::topology<ball, interaction, fem_ds, storage::pattern::empty_item,
                              rt_rt_interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct simulation : simul::time_stepping<td, osi, osnspb> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 1>>;
};
struct make
    : storage::make_runtime_database<
          env, fem_ds, simulation,
          storage::with_properties<
              storage::wrapped<config::interaction, storage::some::unbounded_collection>,
              storage::wrapped<config::relation, storage::some::unbounded_collection>,
              storage::wrapped<config::rt_rt_interaction, storage::some::unbounded_collection>,
              storage::wrapped<config::rt_relation, storage::some::unbounded_collection>,
              storage::unbounded<storage::attr_t<fem_ds, "q">>,
              storage::unbounded<storage::attr_t<fem_ds, "velocity">>,
              storage::unbounded<storage::attr_t<fem_ds, "fext">>,
              storage::sparse<storage::attr_t<fem_ds, "mass_matrix">>,
              storage::sparse<storage::attr_t<fem_ds, "k_matrix">>,
              storage::dynamic_storage<
                  storage::attached<config::ball, storage::pattern::symbol<"bc_velocities_0">,
                                    storage::some::unbounded_vector<storage::some::indice>>>,
              storage::dynamic_storage<storage::attached<
                  config::fem_ds, storage::pattern::symbol<"bc_velocities_0">,
                  storage::some::unbounded_collection<storage::some::indice>>>>> {};

}  // namespace siconos::vkernel::config

namespace {

class FemT3Test : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(FemT3Test);
  CPPUNIT_TEST(testFemT3);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testFemT3() {
    double Ly = 1.0;

    auto gmsh_filename = "square_200.msh";

    // Applied forces
    siconos::algebra::SiconosVector nodal_forces{2};
    nodal_forces << 0., -1e7;
    // Boundary Conditions
    std::vector<int> node_dof_index(2);
    node_dof_index[0] = 0;
    node_dof_index[1] = 1;

    siconos::mechanics::fem::Tags tags;
    tags[siconos::mechanics::fem::MeshTags::bulk_material] = 1;
    tags[siconos::mechanics::fem::MeshTags::boundary_conditions] = 2;
    tags[siconos::mechanics::fem::MeshTags::applied_forces] = 3;

    siconos::mechanics::fem::Material mat{7800, 210e9, 1. / 3};

    auto FEsolid = siconos::mechanics::fem::build_dynamicalsystem_from_gmsh(
        gmsh_filename, tags, mat, nodal_forces, node_dof_index);

    double t0 = 0;     // initial computation time
    double T = 1e-02;  // final computation time
    auto femodel = FEsolid->FEModel();
    double e = 0.0;
    const double gap = Ly * 5e-4;
    siconos::algebra::SiconosVector initial_gap{1};
    initial_gap << gap;
    std::vector<siconos::algebra::SiconosDenseMatrix> Hv;

    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;
    using storage::handle;

    auto data = siconos::vkernel::config::make();
    handle fe_solid = storage::add<config::fem_ds>(data);
    handle nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = e;

    handle lcp = storage::add<config::lcp>(data);
    lcp.create(SICONOS_LCP_LEMKE);

    double h = 1e-05;    // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator

    handle simul = storage::add<config::simulation>(data);
    simul.one_step_integrator().theta() = theta;
    simul.one_step_integrator().constraint_activation_threshold() = 0.;
    simul.time_discretization().t0() = t0;
    simul.time_discretization().h() = h;
    simul.time_discretization().tmax() = T;

    handle osnspb = simul.one_step_nonsmooth_problem();
    osnspb.problem() = lcp;

    handle so = storage::add<siconos::vkernel::simul::solver_options>(data);
    so.create(SICONOS_LCP_LEMKE);
    osnspb.options() = so;

    // the vertical dofs of the bottom nodes, each one in contact with the
    // ground through y = q_y + gap >= 0
    std::vector<std::size_t> contact_dofs;
    for (auto node : femodel->nodes()) {
      if (fabs(node->y()) <= 1e-16 and fabs(node->x()) >= 1e-16) {
        auto idx_y = node->global_dof_index()[1];
        contact_dofs.push_back(idx_y);
        Hv.emplace_back(1, FEsolid->dimension());
        Hv.back().setZero();
        Hv.back()(0, idx_y) = 1.0;

        handle rt_rel = storage::add<config::rt_relation>(data);
        rt_rel.h_matrix().resize(1, FEsolid->dimension());
        rt_rel.h_matrix() = Hv.back();
        rt_rel.b() = initial_gap;

        handle inter = simul.topology().link(fe_solid);
        inter.nslaw() = nslaw;
        simul.topology().bind(inter, rt_rel);
      }
    }
    CPPUNIT_ASSERT_EQUAL_MESSAGE("contact nodes", std::size_t{8}, contact_dofs.size());

    fe_solid.mass_matrix() = FEsolid->mass();
    fe_solid.k_matrix() = FEsolid->stiffnessMatrix();
    fe_solid.fext() = FEsolid->fext();
    fe_solid.velocity() = *FEsolid->velocity();
    fe_solid.q() = *FEsolid->q();

    auto indices = FEsolid->boundaryConditions()->velocityIndices();
    static_assert(storage::is_dynamic_storage_v<config::fem_ds, std::decay_t<decltype(data)>,
                                                "bc_velocities_0">);
    auto& bc_set = storage::prop_values<config::fem_ds, "bc_velocities_0">(data, 0);
    auto& bc_vel = bc_set.emplace(fe_solid.index().value());
    bc_vel.resize(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
      bc_vel[i] = indices[i];
    }

    simul.initialize();

    auto ds_index = storage::prop<"index", config::topo>(fe_solid);
    std::ofstream result_file("fem_t3.dat");
    double min_y = std::numeric_limits<double>::infinity();
    bool finite = true;

    // the columns of T3.cpp: time, q and v of the last dof and of dof 169,
    // the impulse of the first contact
    auto output = [&](double p) {
      auto step = simul.current_step();
      auto q = storage::attr<"q">(fe_solid, step);
      auto v = storage::attr<"velocity">(fe_solid, step);
      auto dim = q.size();
      const double row[] = {
          step * simul.time_step(), q(dim - 1), v(dim - 1), p, q(169), v(169)};
      for (double x : row) finite = finite && std::isfinite(x);
      for (auto dof : contact_dofs) min_y = std::min(min_y, q(dof) + gap);
      std::print(result_file, "{:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n", row[0],
                 row[1], row[2], row[3], row[4], row[5]);
    };

    output(0.);
    while (simul.has_next_event()) {
      unsigned int ninvds = simul.compute_one_step();
      double p = 0.0;
      if (ninvds > 0) {
        handle osi = simul.one_step_integrator();
        p = get_vector(osi.p0_vector_assembled(), ds_index)(0);
      }
      output(p);
    }
    result_file.close();

    const auto dat = siconos::vkernel::test::read_ref("fem_t3.dat");
    CPPUNIT_ASSERT_EQUAL_MESSAGE("1000 steps", Eigen::Index{1001}, dat.rows());
    CPPUNIT_ASSERT_MESSAGE("a finite trajectory", finite);

    // the AOS kernel: the same model, another code, the whole run
    const auto aos = siconos::vkernel::test::read_ref("T3_square_200.ref");
    const auto whole = siconos::vkernel::test::relative_error(dat, aos);
    std::print("{} rows against the AOS reference: relative error {:.3e}\n", dat.rows(),
               whole.relative);
    CPPUNIT_ASSERT_MESSAGE("the AOS reference on the whole run",
                           whole.relative <= aos_tolerance);

    // the ground holds the bottom nodes, up to the penetration of a
    // velocity-level scheme
    std::print("smallest contact gap {:.3e} (gap {:.1e})\n", min_y, gap);
    CPPUNIT_ASSERT_MESSAGE("no penetration beyond the tolerance", min_y >= -penetration);

    CPPUNIT_ASSERT_REF(dat, "fem_t3.ref", 1e-8);
  }

 private:
  static constexpr double aos_tolerance = 1e-12;
  static constexpr double penetration = 1e-4;
};

CPPUNIT_TEST_SUITE_REGISTRATION(FemT3Test);

}  // namespace
