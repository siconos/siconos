#include "siconos/utils/check.hpp"
#include "siconos/config/environment.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/topology.hpp"
#include "siconos/simul/interaction.hpp"
#include "siconos/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/simul/one_step_integrator.hpp"
#include "siconos/simul/time_discretization.hpp"
#include "siconos/simul/time_stepping.hpp"
#include <charconv>
#include <tuple>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <functional>
#include <typeinfo>
#if defined(__clang__)
#include <boost/hana/experimental/type_name.hpp>
#endif
using fmt::print;

using namespace siconos;
using env = siconos::standard_environment;


using namespace boost::hana::literals;
int main()
{
  using ball = lagrangian_ds;
  using lcp = numerics::nonsmooth_problem<LinearComplementarityProblem>;
  using osnspb = numerics::one_step_nonsmooth_problem<lcp>;
  using relation = lagrangian_r;
  using nslaw = nonsmooth_law::newton_impact;
  using interaction = interaction<nslaw, relation, 1>;
  using osi = one_step_integrator<ball, interaction>::moreau_jean;
  using td = time_discretization<>;
  using topo = topology<ball, interaction>;
  using simulation = time_stepping<td, osi, osnspb, topo>;
  using siconos::get;

  struct zz {};

  auto data = siconos::make_storage<env, simulation,
                                    wrap<some::unbounded_collection, relation>,
                                    wrap<some::unbounded_collection, ball>,
                                    wrap<some::unbounded_collection, interaction>,
                                    with_properties<
                                      keep<ball::q, 10>,
                                      diagonal<ball::mass_matrix>,
                                      unbounded_diagonal<osi::mass_matrix_assembled>>>();

  auto the_ball = add<ball>(data);

  the_ball.fext() = { 10., 0., 0.};
  the_ball.q() = { 0., 0., 1.};

   auto& v0 = ball::velocity::at(the_ball);

   v0 = { 1., 2., 3.};

   auto ds1 = add<ball>(data);
   ball::fext::at(the_ball) = { 10., 0., 0.};
   ball::q::at(ds1) = { 1.,1., 1.};
   ball::fext::at(the_ball) = { 10., 0., 0.};

   auto ds2 = add<ball>(data);
   ball::fext::at(the_ball) = { 10., 0., 0.};
   ball::q::at(ds2) = { 9.,9., 9.};
   ball::fext::at(the_ball) = { 10., 0., 0.};
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

//   remove_vertex_item(ds1, 0);
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

   //rint("-----{}------", get<ball::mass_matrix>(ds1, data));
   //print("{}", get<ball::mass_matrix>(the_ball, data));

//   auto m = get<ball::mass_matrix>(the_ball, data);

//   print("---\n");

//   print("{}", m);

//   print("---\n");

   ball::fext::at(the_ball) = { 10., 0., 0.};
//   print("{}\n", data._collections);

   ball::fext::at(the_ball) = { 2.,1.,0.};
//   print("{}\n", data(get<ball::fext>)(the_ball));

   auto inter1 = add<interaction>(data);
   auto nslaw2 = add<nslaw>(data);

   auto dsimul = add<simulation>(data);

   // siconos::set<interaction::nonsmooth_law>(inter1);
   interaction::nonsmooth_law::at(inter1) = nslaw2;
   get<interaction::nonsmooth_law>(inter1, data) = nslaw2;

   auto& e = siconos::get_memory<attr_t<nslaw, "e">>(data);
   e[0][nslaw2.get()] = 0.7;

   nslaw2.e() = 0.6;

   print("{}\n", siconos::get_memory<attr_t<nslaw, "e">>(data));

   auto inter2 = add<interaction>(data);
   auto nslaw3 = add<nslaw>(data);

   interaction::nonsmooth_law::at(inter2) = nslaw3;
   attr_t<nslaw, "e">::at(interaction::nonsmooth_law::at(inter2), data) = 0.3;

   auto htopo = dsimul.topology();
   auto hosi = dsimul.one_step_integrator();
   auto ball1 = add<ball>(data);
   auto ball2 = add<ball>(data);
   auto ball3 = add<ball>(data);
   print("e={}\n", attr_t<nslaw, "e">::at(nslaw3));

   auto intera = htopo.link(ball1);
   auto interb = htopo.link(ball1, ball2);
   auto interc = htopo.link(ball1, ball3);
   auto interd = htopo.link(ball3);

   print("intera.nds={}\n", intera.property(symbol<"nds">{}));
   print("interb.nds={}\n", interb.property(symbol<"nds">{}));
   print("interc.nds={}\n", interc.property(symbol<"nds">{}));
   print("interd.nds={}\n", interd.property(symbol<"nds">{}));

   handle(htopo.dynamical_system_graphs()[0].bundle(ball1.property(symbol<"vd">{})), data).property(symbol<"index">{}) +=10;;

   ball2.property(symbol<"index">{}) +=1;
   print("ball1.index = {}\n", ball1.property(symbol<"index">{}));
   print("ball2.index = {}\n", ball2.property(symbol<"index">{}));

   print("htopo.index : {}\n", mp::get<attached_storage<ball, symbol<"index">, some::indice>>(data));
   print("make_index\n");
   htopo.make_index();
   print("htopo.index : {}\n", mp::get<attached_storage<ball, symbol<"index">, some::indice>>(data));

   hosi.assemble_h_matrix_for_involved_ds(0);
   hosi.assemble_mass_matrix_for_involved_ds(0);

   hosi.compute_w_matrix(0);
   hosi.compute_q_vector_assembled(0);

   auto so = add<numerics::solver_options>(data);
   so.create();

   auto xso = so.instance();
   print("solver options iSize = {}\n", xso->iSize);
   print("SICONOS_IPARAM_MAX_ITER = {}\n", SICONOS_IPARAM_MAX_ITER);
   print("iparam[SICONOS_IPARAM_MAX_ITER] = {}\n", xso->iparam[SICONOS_IPARAM_MAX_ITER]);

   so.instance()->iparam[SICONOS_IPARAM_MAX_ITER] = 100;

   auto lcp_1 = add<lcp>(data);
   lcp_1.create();

   auto osnspb_1 = dsimul.one_step_nonsmooth_problem();
   osnspb_1.problem() = lcp_1;
   osnspb_1.options() = so;

   print("memory_size={}\n", (memory_size<ball::q, decltype(all_properties_as<property::keep>(data))>));

   hosi.compute_q_vector_assembled(0);
   dsimul.solve_nonsmooth_problem<LinearComplementarityProblem>();

   dsimul.compute_input();

   hosi.update_velocity_for_involved_ds();
   hosi.update_all_velocities(0);
   hosi.update_positions(0, 0.01);
//   auto& q = siconos::get_memory<ball::q>(data);
//   print("q={}\n", q);
//   print("v={}\n", siconos::get_memory<ball::velocity>(data));

   auto hsim = make_full_handle<simulation>(0, data);

}
