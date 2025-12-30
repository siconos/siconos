#include "siconos/utils/check.hpp"
#include "siconos/utils/environment.hpp"
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

  auto data0 = siconos::make_storage<env, nslaw>();
  auto nslaw0 = add<nslaw>(data0);

  get<attr_t<nslaw, "e">>(nslaw0, data0) = 0.6;
  get<attr_t<nslaw, "e">>(nslaw0, data0) = 0.7;
  assert(get<attr_t<nslaw, "e">>(nslaw0, data0) == 0.7);

  auto x = get<attr_t<nslaw, "e">>(nslaw0, data0);
  x = 0.4;

  attr_t<nslaw, "e">::at(nslaw0) = 0.8;

  assert (get<attr_t<nslaw, "e">>(nslaw0, data0) == 0.8);

  auto data1 = siconos::make_storage<env, wrap<some::unbounded_collection, nslaw>>();

  auto nslaw1 = add<nslaw>(data1);

  attr_t<nslaw, "e">::at(nslaw1) = 0.9;

  auto data2 = siconos::make_storage<env,
                                     wrap<some::unbounded_collection, nslaw>,
                                     with_properties<attached_storage<nslaw, zz, some::scalar>>>();

  auto nslaw2b = add<nslaw>(data2);

  nslaw2b.property(zz{}) = 2;

  auto data = siconos::make_storage<env, simulation,
                                    wrap<some::unbounded_collection, relation>,
                                    wrap<some::unbounded_collection, ball>,
                                    wrap<some::unbounded_collection, interaction>,
                                    with_properties<
                                      keep<ball::q, 10>,
                                      diagonal<ball::mass_matrix>,
                                      unbounded_diagonal<osi::mass_matrix_assembled>>>();
//   add<simulation>(data);
  print("v={}\n", siconos::get_memory<ball::velocity>(data));
  print("---\n");
  mp::for_each(data,
    [](auto& pair)
    {
      using value_t = std::decay_t<decltype(mp::second(pair))>;

      if constexpr (match::size<value_t>)
      {
        print("->{}\n", mp::second(pair).size());
      }
    });

   print("---\n");


   auto ds0 = add<ball>(data);

//   mp::type_trace<decltype(ball::fext::at(ds0))>();
   ball::fext::at(ds0) = { 10., 0., 0.};
   ball::q::at(ds0) = { 0., 0., 1.};

   auto& v0 = ball::velocity::at(ds0);

   v0 = { 1., 2., 3.};
   print("v={}\n", siconos::get_memory<ball::velocity>(data));
   auto& velocityp = get<ball::velocity>(ds0, data);
   print("--->{}\n", velocityp);

   print("---\n");
   mp::for_each(data,
    [](auto& pair)
    {
      using value_t = std::decay_t<decltype(mp::second(pair))>;

      if constexpr (match::size<value_t>)
      {
        print("->{}\n", mp::second(pair).size());
      }
    });
   print("---\n");

   auto ds1 = add<ball>(data);
   ball::fext::at(ds0) = { 10., 0., 0.};
   ball::q::at(ds1) = { 1.,1., 1.};
   ball::fext::at(ds0) = { 10., 0., 0.};

   auto ds2 = add<ball>(data);
   ball::fext::at(ds0) = { 10., 0., 0.};
   ball::q::at(ds2) = { 9.,9., 9.};
   ball::fext::at(ds0) = { 10., 0., 0.};
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

//   remove_vertex_item(ds1, 0);
//   print("---\n");
//   for_each([](auto& a) { print("{0}\n", a); }, data._collections);

   //rint("-----{}------", get<ball::mass_matrix>(ds1, data));
   //print("{}", get<ball::mass_matrix>(ds0, data));

//   auto m = get<ball::mass_matrix>(ds0, data);

//   print("---\n");

//   print("{}", m);

//   print("---\n");

   ball::fext::at(ds0) = { 10., 0., 0.};
//   print("{}\n", data._collections);

   ball::fext::at(ds0) = { 2.,1.,0.};
//   print("{}\n", data(get<ball::fext>)(ds0));

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

   auto ustore1 = unit_storage<env, some::scalar, some::vector<some::scalar, some::indice_value<3>>, some::graph<some::indice, some::indice>>::type {};

//   mp::find(ustore1, boost::hana::type_c<some::scalar>) = 1.0;
//   mp::find(ustore1, boost::hana::type_c<some::vector<3>>).value() = { 1.0, 2.0, 3.0 };

   mp::get<some::scalar>(ustore1) = 2.0;
   mp::get<some::vector<some::scalar, some::indice_value<3>>>(ustore1) =  {1.0, 2.0, 3.0};

   print("ustore1 scalar={}\n", mp::get<some::scalar>(ustore1));

   print("ustore1 vector={}\n", mp::get<some::vector<some::scalar, some::indice_value<3>>>(ustore1));

   auto ustore2 = typename item_storage<env, simulation, ball>::type {};

   mp::get<ball::q>(ustore2) = { 2.0, 3.0, 4.0 };
   print("ustore2 ball::q ={}\n", mp::get<ball::q>(ustore2));

   //mp::transform(iget<ball>(ustore2), [&ustore2]<typename A>(A){ return iget<A>(ustore2); });

   auto ustore5 = make_storage<env, simulation,
                               wrap<some::unbounded_collection, ball>,
                               wrap<some::unbounded_collection, interaction>,
                               with_properties<diagonal<ball::mass_matrix>,
                                               unbounded_diagonal<osi::mass_matrix_assembled>,
                                               attached_storage<ball, zz , some::scalar>>>();


   //mp::get<ball::q>(ustore5)[0].push_back({7.,8.,9.});
   print("ustore5 ball::q ={}\n", mp::get<ball::q>(ustore5));

   auto some_iball=add<ball>(ustore5);
   get<ball::velocity>(some_iball, ustore5) = {3., 3., 3.};
   print("ustore5 ball::velocity = {}\n", get<ball::velocity>(some_iball, ustore5));

   // ok, compilation failure:
   // get<attr_t<nslaw, "e">>(0, some_iball, ustore5);

   // TD: remove, nsds::link


   get<ball::fext>(some_iball, ustore5) = { 0., 0., 9.81 };

   auto the_td = add<td>(ustore5);
   get<td::step>(the_td, ustore5) = 0;
   get<td::h>(the_td, ustore5) = 0.01;

   auto m = 1.0;
   auto R = 1.0;
//   print("mass matrix: {}", get<ball::mass_matrix>(some_iball, ustore5));
   auto ma = some_iball.mass_matrix();//get<ball::mass_matrix>(some_iball, ustore5);
//   ma.resize(3,3);
   ma.diagonal()[0] = m;
   ma.diagonal()[1] = m;
   ma.diagonal()[2] = 2.5*m*R*R;
//   print("mass matrix: {}", get<ball::mass_matrix>(some_iball, ustore5));

   auto simul = make_full_handle<simulation>(0, ustore5);

   simul.update_indexsets(0);
   simul.compute_one_step();

   print("current step : {}\n", simul.current_step());
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));
   simul.compute_one_step();
   print("current step : {}\n", simul.current_step());
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));
   simul.compute_one_step();
   print("current step : {}\n", simul.current_step());
   print("v={}\n", siconos::get_memory<ball::velocity>(ustore5));

   auto a = get<ball::velocity>(simul.current_step(), some_iball, ustore5);
   print("ustore5 ball::velocity = {}\n", a);

   auto& dsg0 = get<topo::dynamical_system_graphs>(htopo, ustore5)[0];

   auto dsgv = dsg0.add_vertex(some_iball);

   auto& ig0 = get<topo::interaction_graphs>(htopo, ustore5)[0];

   auto inter = add<interaction>(ustore5);

   auto [new_edge, ig_new_ve] = dsg0.add_edge(dsgv, dsgv, inter, ig0);
   auto rel = add<relation>(ustore5);

   auto yyv = mp::get<interaction::h_matrix>(ustore5);
   interaction::h_matrix::at(inter)(0, 0) = 1;

   interaction::relation::at(inter) = rel;
   auto& v = interaction::y::at(inter)[0];
   v = { 1 };

//   rel.compute_output(0., some_iball, inter, 1_c);

//   rel.compute_input(0., some_iball, inter, 1_c);
   auto hsim = make_full_handle<simulation>(0, data);
//   hsim.compute_output(1_c);

   mp::get<attached_storage<ball, zz, some::scalar>>(ustore5)[0][0] = 1.0;
   get<zz>(0, some_iball, ustore5) = 1.0;
   some_iball.property(zz{}) = 1.0;
//   for_each(ustore5,
   //           []<typename Key, typename Value>(Key&& key, Value&& value)
   //        {
//              print("{}\n", boost::hana::experimental::type_name<Key>().c_str());
   //       });


#ifdef __clang__
//   print("-->{}\n", boost::hana::experimental::type_name<mp::pair<ball, zz>>().c_str());
#endif
   auto dd = make_storage<env, wrap<some::unbounded_collection, ball>,
                          with_properties<diagonal<ball::mass_matrix>, keep<ball::q, 10>,
                                          attached_storage<ball, decltype("z"_s), some::scalar>,
                                          attached_storage<ball, symbol<"x">, some::scalar>>>();

  auto xball1 = add<ball>(dd);
  auto xball2 = add<ball>(dd);
  auto xball3 = add<ball>(dd);

  print("{},{}\n", xball1.get(), xball2.get());
  xball1.velocity() = {1, 1, 1};
  xball2.velocity() = {2, 2, 2};
  xball3.velocity() = {3, 3, 3};
  print("{}\n", mp::get<ball::velocity>(dd));
  siconos::remove(xball1, dd);
  print("{},{},{}\n", xball1.get(), xball2.get(), xball3.get());
  print("{}\n", mp::get<ball::velocity>(dd));

  static_assert(match::tag<attached_storage<ball, zz, some::scalar>, zz>);
  auto xball4 = add<ball>(dd);
  xball4.property(symbol<"x">{}) = 3.0;
  xball4.property("z"_s) = 4.0;
  xball2.property(symbol<"x">{}) = 1.0;
  xball2.property("z"_s) = 2.0;



  print("{},{},{},{}\n", xball2.property(symbol<"x">{}), xball2.property("z"_s), xball4.property(symbol<"x">{}), xball4.property("z"_s));

   for_each(dd, []<typename Key, typename Value>(Key k, Value v)
            {
              print("Key: [{}]\n", mp::type_name<Key>());
              print("Value: [{}]\n\n", mp::type_name<Value>());
            });
//   assert((mp::get<attached_storage<some::unbounded_collection<ball>, internal_index, some::indice>>(dd)[0][1] == 1));

}
