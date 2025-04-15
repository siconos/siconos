#include <concepts>

#include "siconos/config/config.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/utils/print.hpp"

namespace siconos::config {
using params = map<assoc<param<"dof">, param_val<3>>>;
}

using namespace siconos;

using env = standard_environment<config::params>;

#include "siconos/collision/point.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/neighborhood.hpp"
#include "siconos/collision/space_filter.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/traits/traits.hpp"
#include "siconos/utils/environment.hpp"

namespace siconos {

namespace some = storage::some;
namespace traits = storage::traits;
namespace ground = storage::ground;
using namespace storage::pattern;

struct aaa {
  int v = 1;
};
struct bbb : item<> {
  struct attr : some::specific<pointer<aaa>>, siconos::access<attr> {};
  using attributes = gather<attr>;
  template <typename H>
  struct interface : siconos::default_interface<H> {
    constexpr int hello() { return 36; };

    auto methods()
    {
      return collect(method(symbol<"hello">{}, &interface<H>::hello));
    };
  };
};

static_assert(std::is_same_v<traits::config<standard_environment<int>>::
                                 convert<some::scalar>::type,
                             double>);
static_assert(std::is_same_v<traits::config<standard_environment<int>>::
                                 convert<pointer<float>>::type,
                             pointer<float>>);

static_assert(
    std::is_same_v<traits::config<standard_environment<int>>::convert<
                       some::specific<pointer<float>>>::type,
                   pointer<float>>);

static_assert(
    match::abstract_matrix<
        attribute<"attr0", some::matrix<some::scalar, some::indice_value<1>,
                                        some::indice_value<1>>>>);

struct item0 : item<> {
  using attributes = gather<
      attribute<"attr0", some::matrix<some::scalar, some::indice_value<1>,
                                      some::indice_value<1>>>>;
};

static_assert(
    boost::hana::is_convertible<
        siconos::storage::pattern::paired<
            siconos::item0,
            siconos::storage::pattern::attribute<
                "attr0", siconos::storage::some::matrix<
                             siconos::storage::some::scalar,
                             siconos::storage::some::indice_value<1>,
                             siconos::storage::some::indice_value<1>>>>,
        siconos::storage::pattern::paired<
            siconos::item0,
            siconos::storage::pattern::attribute<
                "attr0", siconos::storage::some::matrix<
                             siconos::storage::some::scalar,
                             siconos::storage::some::indice_value<1>,
                             siconos::storage::some::indice_value<1>>>>>::
        value);

static_assert(
    std::is_same_v<
        decltype(storage::pattern::attributes(item0{})),
        gather<storage::pattern::paired<
            item0, attribute<"attr0",
                             some::matrix<some::scalar, some::indice_value<1>,
                                          some::indice_value<1>>>>>>);

static_assert(
    std::is_same_v<
        typename std::decay_t<decltype(ground::get<storage::info>(
            storage::make<standard_environment<int>, item0,
                          storage::with_properties<
                              storage::diagonal<attr_t<item0, "attr0">>>>()
                .store()))>::all_properties_t,
        gather<siconos::storage::diagonal<siconos::storage::pattern::paired<
            siconos::item0,
            siconos::storage::pattern::attribute<
                "attr0", siconos::storage::some::matrix<
                             siconos::storage::some::scalar,
                             siconos::storage::some::indice_value<1>,
                             siconos::storage::some::indice_value<1>>>>>>>);

// static_assert(ground::filter(
//     ground::transform(
//         typename std::decay_t<decltype(ground::get<storage::info>(
//             storage::make<standard_environment<int>, item0,
//                           storage::with_properties<storage::diagonal<attr_t<
//                               item0, "attr0">>>>()))>::all_properties_t{},
//         []<typename P>(P) { return typename P::type{}; }),
//     ground::is_parent<siconos::storage::pattern::paired<
//         siconos::item0,
//         siconos::storage::pattern::attribute<
//             "attr0", siconos::storage::some::matrix<
//                          siconos::storage::some::scalar,
//                          siconos::storage::some::indice_value<1>,
//                          siconos::storage::some::indice_value<1>>>>>));

// ground::filter(typename info_t::all_properties_t{},
//                      ground::is_a_model<[]<typename T>() constexpr {
//                        return match::attached_storage<T, Item>;
//                      }>);

static_assert(
    std::is_same_v<
        decltype(ground::filter(
            typename std::decay_t<decltype(ground::get<storage::info>(
                storage::make<standard_environment<int>, item0,
                              storage::with_properties<storage::attached<
                                  item0, symbol<"zz">, some::scalar>>>()
                    .store()))>::all_properties_t{},
            ground::is_a_model<[]<typename T>() constexpr {
              return match::attached_storage<T, item0>;
            }>)),
        ground::tuple<siconos::storage::attached<
            siconos::item0, siconos::storage::pattern::symbol<"zz">,
            siconos::storage::some::scalar>>>);

static_assert(match::diagonal_matrix<
              decltype(storage::attr<"attr0">(storage::add<item0>(
                  storage::make<standard_environment<int>, item0,
                                storage::with_properties<storage::diagonal<
                                    attr_t<item0, "attr0">>>>())))>);

static_assert(
    match::mat<
        std::decay_t<decltype(storage::attr<"attr0">(storage::add<item0>(
            storage::make<standard_environment<int>, item0,
                          storage::with_properties<storage::assembled_matrix<
                              attr_t<item0, "attr0">>>>())))>>);

/* alias : same type, item1 has a diagonal matrix also */
using item1 = item0;

static_assert(match::diagonal_matrix<
              decltype(storage::attr<"attr0">(storage::add<item1>(
                  storage::make<standard_environment<int>, item0, item1,
                                storage::with_properties<storage::diagonal<
                                    attr_t<item0, "attr0">>>>())))>);

/* new type, item2 does not have a diagonal matrix */
struct item2 : item0 {};
static_assert(!match::diagonal_matrix<
              decltype(storage::attr<"attr0">(storage::add<item2>(
                  storage::make<standard_environment<int>, item0, item2,
                                storage::with_properties<storage::diagonal<
                                    attr_t<item0, "attr0">>>>())))>);

static_assert(
    match::matrix<decltype(storage::attr<"attr0">(storage::add<item0>(
        storage::make<standard_environment<int>,
                      wrap<some::unbounded_collection, item0>>())))>);

static_assert(
    match::push_back<
        std::decay_t<decltype(storage::attr_values<item0, "attr0">(
            storage::add<item0>(
                storage::make<standard_environment<int>,
                              wrap<some::unbounded_collection, item0>>())
                .data(),
            0))>>);

static_assert(match::diagonal_matrix<
              decltype(storage::attr<"attr0">(storage::add<item0>(
                  storage::make<standard_environment<int>,
                                wrap<some::unbounded_collection, item0>,
                                storage::with_properties<storage::diagonal<
                                    attr_t<item0, "attr0">>>>())))>);

static_assert(
    std::is_same_v<
        std::decay_t<decltype((storage::add<model::newton_impact>(
                                   storage::make<standard_environment<int>,
                                                 model::newton_impact>()))
                                  .e())>,
        typename standard_environment<int>::scalar>);

static_assert(
    std::is_same_v<
        std::decay_t<
            decltype(storage::add<collision::shape::disk>(
                         storage::make<
                             standard_environment<int>,
                             collision::shape::disk,
                             storage::with_properties<storage::bind<
                                 collision::shape::disk, "disk_shape">>>())
                         .radius())>,
        typename standard_environment<int>::scalar>);

//   error: static assertion expression is not an integral constant expression
// static_assert(
//     storage::has_property<collision::shape::disk, storage::property::bind>(
//         storage::make<standard_environment<int>, collision::shape::disk,
//                       storage::with_properties<storage::bind<
//                           collision::shape::disk, "disk_shape">>>()));

}  // namespace siconos

using namespace boost::hana::literals;

using ball = model::lagrangian_ds;
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = simul::one_step_nonsmooth_problem<lcp>;
using nslaw = model::newton_impact;
using relation = model::lagrangian_r<nslaw::size>;
using interaction = simul::interaction<nslaw, relation>;
using topo = simul::topology<ball, interaction>;
using osi = simul::one_step_integrator<topo>::moreau_jean;
using td = simul::time_discretization<>;
using disk = collision::shape::disk;
using simulation = simul::time_stepping<td, osi, osnspb>;

using disk_shape = collision::shape::disk;

// static_assert(
//   storage::attr<"nslaws">(
//     storage::add<inter_manager>(
//         storage::make <
//         standard_environment<config::map<config::iparam<"ncgroups",
// 1>>>,
//         inter_manager, nslaw>())
//         .insert_nonsmooth_law(
//             storage::add<nslaw>(
//                 storage::make <
// standard_environment<config::map<config::iparam<"ncgroups",
//                 1>>>, inter_manager, nslaw>()),
//             0, 0))(0, 0) == storage::add<nslaw>(
//                 storage::make <
// standard_environment<config::map<config::iparam<"ncgroups",
//                 1>>>, inter_manager, nslaw>()));

static_assert(
    match::diagonal_matrix<
        decltype(storage::attr<"mass_matrix">(storage::add<ball>(
            storage::make<
                standard_environment<config::map<
                    config::iparam<"dof", 3>, config::iparam<"ncgroups", 1>>>,
                simulation, ball, relation, interaction,
                storage::with_properties<
                    storage::diagonal<attr_t<ball, "mass_matrix">>>>())))>);

static_assert(
    match::index<std::decay_t<decltype(storage::prop<
                                       "shape">(storage::add<ball>(
        storage::make<
            standard_environment<config::map<config::iparam<"dof", 3>>>, ball,
            storage::with_properties<storage::attached<
                ball, symbol<"shape">, some::item_ref<disk_shape>>>>())))>>);

template <typename T>
struct is_polymorhic : std::integral_constant<bool, []() {
  return match::polymorphic_type<T>;
}()> {};

static_assert(match::relation1<storage::handle<
                  relation, int, decltype(storage::make<env, relation>())>>);

//  {
static_assert(std::is_same_v<decltype(all_items(nslaw{})),
                             gather<siconos::model::newton_impact>>);

static_assert(
    std::derived_from<std::decay_t<decltype(ground::filter(
                          typename interaction::attributes{},
                          ground::compose(ground::trait<is_polymorhic>,
                                          ground::typeid_))[0_c])>,
                      some::polymorphic_attribute<some::item_ref<relation>>>);

// ground::type_trace<decltype(all_items(interaction{}))>();
static_assert(std::is_same_v<decltype(all_items(interaction{})),
                             gather<interaction, nslaw, relation>>);

// static_assert(must::contains<osnspb,
// decltype(all_items(simulation{}))>);

static_assert(match::item<ball>);
static_assert(match::attribute<attr_t<nslaw, "e">>);

static_assert(match::attribute_of<attr_t<nslaw, "e">, nslaw>);

static_assert(match::attribute_of<attr_t<ball, "velocity">, ball>);
static_assert(match::attribute_of<attr_t<td, "step">, td>);

// static_assert(std::is_same_v<td,
// decltype(item_attribute<td::step>(all_items(simulation{})))>);

static_assert(
    std::is_same_v<typename siconos::traits::config<env>::convert<attr_t<
                       siconos::simul::time_discretization<>, "step">>::type,
                   typename env::indice>);

static_assert(
    std::is_same_v<
        typename siconos::traits::config<env>::convert<some::scalar>::type,
        typename env::scalar>);

struct zz {};

static_assert(
    match::attached_storage<storage::attached<ball, zz, some::scalar>, ball>);
static_assert(
    match::attached_storage<storage::attached<ball, zz, some::scalar>,
                            wrap<some::unbounded_collection, ball>>);
static_assert(match::attached_storage<
              storage::attached<ball, symbol<"Z">, some::scalar>,
              wrap<some::unbounded_collection, ball>>);
static_assert(!match::attached_storage<
              storage::attached<ball, zz, some::scalar>, nslaw>);

static_assert(
    match::unbounded_storage<some::unbounded_collection<some::scalar>>);
static_assert(match::bounded_storage<
              some::bounded_collection<some::scalar, some::indice_value<1>>>);
static_assert(!match::unbounded_storage<
              some::bounded_collection<some::scalar, some::indice_value<1>>>);
static_assert(
    !match::bounded_storage<some::unbounded_collection<some::scalar>>);
static_assert(
    match::unbounded_storage<some::unbounded_diagonal_matrix<some::scalar>>);

static_assert(match::wrap<wrap<some::unbounded_diagonal_matrix, ball>>);

static_assert(match::bounded_storage<
              wrap<some::bounded_collection, ball,
                   some::indice_value<1>>::template wrapper<some::scalar>>);

static_assert(traits::translatable<int, env>);

static_assert(std::is_same_v<traits::config<env>::convert<int>::type, int>);

static_assert(
    traits::translatable<siconos::some::unbounded_collection<char[3]>, env>);

static_assert(
    ground::filter(gather<char, int, double>{},
                   ground::compose(ground::trait<std::is_floating_point>,
                                   ground::typeid_)) == gather<double>{});

static_assert(std::is_same_v<
              decltype(ground::filter(
                  std::tuple<attr_t<nslaw, "e">, attr_t<ball, "velocity">>{},
                  ground::derive_from<some::scalar>)),
              std::tuple<attr_t<nslaw, "e">>>);

static_assert(
    std::is_same_v<
        decltype(attributes(interaction{})),
        gather<attr_t<interaction, "relation">, attr_t<interaction, "nslaw">,
               attr_t<interaction, "h_matrix1">,
               attr_t<interaction, "h_matrix2">,
               attr_t<interaction, "lambda">, attr_t<interaction, "y">,
               attr_t<interaction, "ydot">>>);

static_assert(
    std::is_same_v<
        decltype(flatten(all_attributes(interaction{}))),
        gather<attr_t<interaction, "relation">, attr_t<interaction, "nslaw">,
               attr_t<interaction, "h_matrix1">,
               attr_t<interaction, "h_matrix2">,
               attr_t<interaction, "lambda">, attr_t<interaction, "y">,
               attr_t<interaction, "ydot">, attr_t<nslaw, "e">,
               attr_t<relation, "b">, attr_t<relation, "h_matrix">>>);

//}  // namespace siconos

namespace siconos::config {
using disk = model::lagrangian_ds;
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = simul::one_step_nonsmooth_problem<lcp>;
using nslaw = model::newton_impact;
using disk_shape = collision::shape::disk;
using diskdisk_r = collision::diskdisk_r;
using diskline_r = collision::diskline_r;
using interaction = simul::interaction<nslaw, diskdisk_r, diskline_r>;
using topo = simul::topology<disk, interaction>;
using osi = simul::one_step_integrator<topo>::moreau_jean;
using td = simul::time_discretization<>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;
using pointd = collision::point<disk>;
using pointl = collision::point<collision::shape::line>;
using neighborhood = collision::neighborhood<pointd, pointl>;
using space_filter = collision::space_filter<topo, neighborhood>;

using params = map<iparam<"dof", 3>>;
}  // namespace siconos::config

// static_assert(
//     typename std::decay_t<decltype(ground::get<storage::info>(
//         storage::make<
//             standard_environment<config::params>, config::simulation,
//             wrap<some::unbounded_collection, config::disk>,
//             config::disk_shape, config::diskdisk_r,
//             wrap<some::unbounded_collection, config::diskline_r>,
//             wrap<some::unbounded_collection, config::pointl>,
//             wrap<some::unbounded_collection, config::pointd>,
//             wrap<some::unbounded_collection, config::interaction>,
//             config::space_filter,
//             storage::with_properties<
//                 storage::attached<
//                     config::disk, storage::pattern::symbol<"shape">,
//                     storage::some::item_ref<config::disk_shape>>,
//                 storage::time_invariant<storage::attr_t<config::disk, "fex"
//                                                                       "t">>,
//                 storage::diagonal<
//                     storage::attr_t<config::disk, "mass_matrix">>,
//                 storage::unbounded_diagonal<storage::attr_t<
//                     config::osi, "mass_matrix_assembled">>>>()))>::
//         all_properties_t{});

int main() {}
