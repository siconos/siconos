// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

// storage::remove keeps the dynamic properties of the items consistent.
//
// Three balls, ball i at q(0) = i.  Balls 0 and 2 get a bc_velocities_0 set
// (a dynamic property attached by the topology).  Removing ball 0 moves ball 2
// into index 0 (swap and pop): its BC set must follow it, and nothing may be
// left at index 2.

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <initializer_list>
#include <map>
#include <print>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "siconos/vkernel.hpp"
#include "siconos/vkernel/config/config_builder.hpp"
#include "siconos/vkernel/storage/permute.hpp"
#include "siconos/vkernel/storage/refs.hpp"

namespace siconos::vkernel::config {
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
                             .build_runtime_database()){};
}  // namespace siconos::vkernel::config

// a per-item storage that cannot grow (a fixed-size array), and an item with
// no dense column at all (a dynamic property only): adding either invents an
// item, so storage::add refuses them
namespace siconos::vkernel::config {
struct fixed_row : storage::pattern::item {
  using attributes =
      storage::pattern::gather<storage::pattern::attribute<"v", storage::some::scalar>>;
};

template <typename... Ts>
struct one_array_helper;
template <typename T>
struct one_array_helper<T> {
  using type = storage::some::array<T, storage::some::indice_value<1>>;
};
template <typename... Ts>
using one_array = typename one_array_helper<Ts...>::type;

struct fixed_rows_t : storage::make_runtime_database<
                          standard_environment, fixed_row,
                          storage::with_properties<storage::wrapped<fixed_row, one_array>>> {};

struct sparse_row : storage::pattern::item {};
struct sparse_rows_t
    : storage::make_runtime_database<
          standard_environment, sparse_row,
          storage::with_properties<storage::dynamic_storage<storage::attached<
              sparse_row, storage::pattern::symbol<"weight">, storage::some::scalar>>>> {};
}  // namespace siconos::vkernel::config

// items without attributes: one counted through its attached storage, one
// without any column
namespace siconos::vkernel::config {
struct marker : storage::pattern::item {};
struct bare : storage::pattern::item {};
struct markers_t
    : storage::make_runtime_database<
          standard_environment, marker, bare,
          storage::with_properties<storage::attached<
              marker, storage::pattern::symbol<"weight">, storage::some::scalar>>> {};
}  // namespace siconos::vkernel::config

// the simulation and its integrator live as long as the store: storage::remove
// refuses them at compile time; a ball may go
namespace match = siconos::vkernel::storage::pattern::match;
static_assert(match::not_removable<siconos::vkernel::config::simulation>);
static_assert(match::not_removable<siconos::vkernel::config::osi>);
static_assert(!match::not_removable<siconos::vkernel::config::ball>);
static_assert(!match::not_removable<siconos::vkernel::config::interaction>);

// a reference handle holds a reference to its index: built from a temporary
// index it would dangle, so that construction must not exist (the same reason
// make_ref_handle(D&, index&&) is deleted)
using ball_ref_handle =
    siconos::vkernel::storage::handle<siconos::vkernel::storage::handle_ref,
                                      siconos::vkernel::config::ball, std::size_t,
                                      siconos::vkernel::config::data_t>;
using ball_ref_index =
    siconos::vkernel::storage::index<siconos::vkernel::config::ball, std::size_t>;
static_assert(!std::constructible_from<ball_ref_handle, siconos::vkernel::config::data_t&,
                                       ball_ref_index&&>);
static_assert(!std::constructible_from<ball_ref_handle, siconos::vkernel::config::data_t,
                                       ball_ref_index&&>);

namespace {

class RegistryRemoveTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(RegistryRemoveTest);
  CPPUNIT_TEST(testRegistryRemove);
  CPPUNIT_TEST(testReadCreatesNoEntry);
  CPPUNIT_TEST(testRemovePastTheEnd);
  CPPUNIT_TEST(testPermuteChecks);
  CPPUNIT_TEST(testPermuteRejectsNonBijection);
  CPPUNIT_TEST(testRefInvalidOfItsOwnType);
  CPPUNIT_TEST(testTableRemapPastTheTable);
  CPPUNIT_TEST(testRemoveRemapsRefs);
  CPPUNIT_TEST(testPermuteRemapsRefs);
  CPPUNIT_TEST(testRemoveManyRemapsRefs);
  CPPUNIT_TEST(testRemapWalker);
  CPPUNIT_TEST(testHandlesWithoutAttributes);
  CPPUNIT_TEST(testFixedSizePerItemStorageIsRefused);
  CPPUNIT_TEST(testSparseOnlyItemIsRefused);
  CPPUNIT_TEST(testConstAccess);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testRegistryRemove() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;

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

    CPPUNIT_ASSERT_MESSAGE("two balls left", std::size(qs) == 2);
    CPPUNIT_ASSERT_MESSAGE("ball 2 moved to index 0", qs[0](0) == 2.);
    CPPUNIT_ASSERT_MESSAGE("ball 1 still at index 1", qs[1](0) == 1.);
    CPPUNIT_ASSERT_MESSAGE("ball 2's BC set followed it to index 0",
                           bcs.contains(0) && is_bc(bcs.at(0), {2, 22}));
    CPPUNIT_ASSERT_MESSAGE("ball 1 still has no BC set", !bcs.contains(1));
    CPPUNIT_ASSERT_MESSAGE("nothing left at index 2", !bcs.contains(2));
  }

  // a read of a dynamic property creates no entry for an item without it: it
  // would change which items the loops over the entries visit, and what is
  // serialized
  void testReadCreatesNoEntry() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;

    auto ball0 = storage::add<config::ball>(data);
    auto ball1 = storage::add<config::ball>(data);

    auto& bc0 = storage::prop<"bc_velocities_0", config::topo>(ball0);
    bc0.resize(1);
    bc0[0] = 7;

    auto& bcs = storage::prop_values<config::ball, "bc_velocities_0", config::topo>(data, 0);
    CPPUNIT_ASSERT_MESSAGE("one BC set", bcs.size() == 1);

    const auto& read1 = storage::read_prop<"bc_velocities_0", config::topo>(ball1);
    CPPUNIT_ASSERT_MESSAGE("ball 1 reads an empty BC set", read1.size() == 0);
    CPPUNIT_ASSERT_MESSAGE("the read of ball 1 created no entry",
                           bcs.size() == 1 && !bcs.contains(1));

    const auto& read0 = storage::read_prop<"bc_velocities_0", config::topo>(ball0);
    CPPUNIT_ASSERT_MESSAGE("ball 0 reads its BC set", read0.size() == 1 && read0[0] == 7);
    CPPUNIT_ASSERT_MESSAGE("still one BC set", bcs.size() == 1);
  }

  // a handle kept across a removal may hold an index past the end (swap and
  // pop): its removal throws instead of moving columns out of bounds
  void testRemovePastTheEnd() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;

    auto ball0 = storage::add<config::ball>(data);
    auto ball1 = storage::add<config::ball>(data);
    auto stale = ball1;
    storage::remove(data, ball1);

    CPPUNIT_ASSERT_THROW(storage::remove(data, stale), std::out_of_range);
    CPPUNIT_ASSERT_MESSAGE("ball 0 still there",
                           std::size(storage::attr_values<config::ball, "q">(data, 0)) == 1);
    storage::remove(data, ball0);
  }

  // storage::permute takes a permutation of the items
  void testPermuteChecks() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;

    auto ball0 = storage::add<config::ball>(data);
    auto ball1 = storage::add<config::ball>(data);
    ball0.q() = {0., 0., 0.};
    ball1.q() = {1., 0., 0.};

    CPPUNIT_ASSERT_THROW(storage::permute<config::ball>(data, {0}), std::invalid_argument);
    CPPUNIT_ASSERT_THROW(storage::permute<config::ball>(data, {0, 5}), std::invalid_argument);

    storage::permute<config::ball>(data, {1, 0});
    auto& qs = storage::attr_values<config::ball, "q">(data, 0);
    CPPUNIT_ASSERT_MESSAGE("swapped", qs[0](0) == 1. && qs[1](0) == 0.);
  }

  // storage::permute takes a permutation of the items, in every build: a
  // duplicated index must be refused, not silently applied
  void testPermuteRejectsNonBijection() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;
    storage::add<config::ball>(data);
    storage::add<config::ball>(data);

    CPPUNIT_ASSERT_THROW(storage::permute<config::ball>(data, {0, 0}), std::invalid_argument);
  }

  // a ref to a removed row becomes the invalid of its own index type, not a
  // truncation of ref_remap::invalid (SIZE_MAX): with an indice narrower than
  // size_t, valid() must still say no item
  void testRefInvalidOfItsOwnType() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;
    using small_ref = storage::index<config::ball, int>;

    const auto kind = storage::kind_id<config::ball>::get();
    const auto r = storage::ref_remap::swap_and_pop(0, 1);  // 0 dies, 1 -> 0

    std::vector<small_ref> refs{small_ref(0), small_ref(1)};
    storage::refs_detail::remap<storage::kind_id>(refs, kind, r);

    CPPUNIT_ASSERT_MESSAGE("the moved ref follows its row", refs[1].value() == 0);
    CPPUNIT_ASSERT_MESSAGE("the removed ref is the index type's invalid",
                           refs[0].value() == small_ref::invalid);
    CPPUNIT_ASSERT_MESSAGE("and valid() says it refers to no item", !refs[0].valid());
  }

  // interactions between balls: inter k links balls a[k] and b[k]
  template <typename Data, typename Balls>
  static auto link(Data& data, Balls& balls,
                   std::initializer_list<std::pair<int, int>> pairs) {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;
    std::vector<decltype(storage::add<config::interaction>(data))> inters;
    for (auto [a, b] : pairs) {
      auto inter = storage::add<config::interaction>(data);
      storage::prop<"ds1", config::topo>(inter) = balls[a];
      storage::prop<"ds2", config::topo>(inter) = balls[b];
      inters.push_back(inter);
    }
    return inters;
  }

  // a table remap covers the rows that existed before the removal: a ref past
  // the table refers to no row, it does not keep its value
  void testTableRemapPastTheTable() {
    namespace storage = siconos::vkernel::storage;

    const std::vector<std::size_t> old_to_new{storage::ref_remap::invalid, 0};
    const auto r = storage::ref_remap::table(old_to_new);

    CPPUNIT_ASSERT_MESSAGE("a row of the table maps to its new row", r(1) == 0);
    CPPUNIT_ASSERT_MESSAGE("a removed row maps to no row",
                           r(0) == storage::ref_remap::invalid);
    CPPUNIT_ASSERT_MESSAGE("a row past the table refers to no row",
                           r(5) == storage::ref_remap::invalid);
  }

  using ends_t = std::pair<std::size_t, std::size_t>;

  template <typename Data>
  static ends_t ends(Data& data, std::size_t k) {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;
    auto& ds1s = storage::prop_values<config::interaction, "ds1", config::topo>(data, 0);
    auto& ds2s = storage::prop_values<config::interaction, "ds2", config::topo>(data, 0);
    return {ds1s[k].value(), ds2s[k].value()};
  }

  static constexpr std::size_t dead = siconos::vkernel::storage::ref_remap::invalid;

  // storage::remove rewrites the refs to the removed item (invalid) and to
  // the moved one (its new index) held in the other columns
  void testRemoveRemapsRefs() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;
    std::vector balls{storage::add<config::ball>(data), storage::add<config::ball>(data),
                      storage::add<config::ball>(data)};
    link(data, balls, {{0, 2}, {1, 2}});

    storage::remove(data, balls[0]);  // ball 2 moves to index 0

    CPPUNIT_ASSERT_MESSAGE("inter 0: ball 0 gone, ball 2 at 0",
                           ends(data, 0) == ends_t(dead, 0));
    CPPUNIT_ASSERT_MESSAGE("inter 1: ball 1 kept, ball 2 at 0", ends(data, 1) == ends_t(1, 0));
  }

  // storage::permute: the refs follow the moved items
  void testPermuteRemapsRefs() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;
    std::vector balls{storage::add<config::ball>(data), storage::add<config::ball>(data),
                      storage::add<config::ball>(data)};
    link(data, balls, {{0, 1}, {1, 2}});

    // new 0 = old 2, new 1 = old 0, new 2 = old 1
    storage::permute<config::ball>(data, {2, 0, 1});

    CPPUNIT_ASSERT_MESSAGE("inter 0: (0, 1) -> (1, 2)", ends(data, 0) == ends_t(1, 2));
    CPPUNIT_ASSERT_MESSAGE("inter 1: (1, 2) -> (2, 0)", ends(data, 1) == ends_t(2, 0));
  }

  // storage::remove_many: several swaps and pops, the refs rewritten once
  void testRemoveManyRemapsRefs() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;
    std::vector balls{storage::add<config::ball>(data), storage::add<config::ball>(data),
                      storage::add<config::ball>(data), storage::add<config::ball>(data),
                      storage::add<config::ball>(data)};
    for (std::size_t k = 0; k < balls.size(); ++k) balls[k].q() = {double(k), 0., 0.};
    link(data, balls, {{0, 1}, {2, 3}, {3, 4}});

    storage::remove_many<config::ball>(data, {1, 3, 1});

    auto& qs = storage::attr_values<config::ball, "q">(data, 0);
    CPPUNIT_ASSERT_MESSAGE("three balls left", std::size(qs) == 3);
    // where the old balls 0, 2, 4 are now
    std::size_t where[5] = {dead, dead, dead, dead, dead};
    for (std::size_t k = 0; k < std::size(qs); ++k) where[std::size_t(qs[k](0))] = k;
    CPPUNIT_ASSERT_MESSAGE("balls 0, 2, 4 kept",
                           where[0] != dead && where[2] != dead && where[4] != dead);
    CPPUNIT_ASSERT_MESSAGE("inter (0, 1)", ends(data, 0) == ends_t(where[0], dead));
    CPPUNIT_ASSERT_MESSAGE("inter (2, 3)", ends(data, 1) == ends_t(where[2], dead));
    CPPUNIT_ASSERT_MESSAGE("inter (3, 4)", ends(data, 2) == ends_t(dead, where[4]));
  }

  // storage::handles counts the items from any of their columns, not only
  // from a first attribute: an item may have none
  void testHandlesWithoutAttributes() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::markers_t data;
    for (int k = 0; k < 3; ++k) storage::add<config::marker>(data);
    CPPUNIT_ASSERT_MESSAGE("three markers",
                           std::ranges::distance(storage::handles<config::marker>(data)) == 3);
    storage::remove(data, *storage::handles<config::marker>(data).begin());
    CPPUNIT_ASSERT_MESSAGE("two markers after a removal",
                           std::ranges::distance(storage::handles<config::marker>(data)) == 2);
    CPPUNIT_ASSERT_MESSAGE("an item without any column has no handle",
                           std::ranges::distance(storage::handles<config::bare>(data)) == 0);
    // a placeholder item (no column at all: an interaction's relation) may be
    // added; its index is the meaningless 0, and removing it is a no-op
    CPPUNIT_ASSERT_NO_THROW(storage::add<config::bare>(data));
    storage::remove(data, *storage::handles<config::marker>(data).begin());
  }

  // a per-item storage that cannot grow would hold every item in slot 0:
  // refuse the add instead of inventing one item
  void testFixedSizePerItemStorageIsRefused() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::fixed_rows_t data;
    CPPUNIT_ASSERT_THROW(storage::add<config::fixed_row>(data), std::logic_error);
    CPPUNIT_ASSERT_MESSAGE(
        "nothing was added",
        std::ranges::distance(storage::handles<config::fixed_row>(data)) == 0);
  }

  // an item whose only column is a dynamic property has no row count: a
  // removal could not follow it, so the add is refused
  void testSparseOnlyItemIsRefused() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::sparse_rows_t data;
    CPPUNIT_ASSERT_THROW(storage::add<config::sparse_row>(data), std::logic_error);
  }

  // a store that is only read can be const: the attribute and property
  // columns are then read through const references, which is what a
  // serializer or a read-only view needs
  void testConstAccess() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;

    config::data_t data;
    auto ball0 = storage::add<config::ball>(data);
    ball0.q() = {1., 2., 3.};
    auto& bc0 = storage::prop<"bc_velocities_0", config::topo>(ball0);
    bc0.resize(1);
    bc0[0] = 7;

    const config::data_t& cdata = data;

    const auto& qs = storage::attr_values<config::ball, "q">(cdata, 0);
    CPPUNIT_ASSERT_MESSAGE("a const store reads q", qs[0](1) == 2.);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(qs[0])>>,
                  "a const store gives const attribute values");

    const auto& bcs =
        storage::prop_values<config::ball, "bc_velocities_0", config::topo>(cdata, 0);
    CPPUNIT_ASSERT_MESSAGE("a const store reads a dynamic property",
                           bcs.contains(0) && bcs.at(0)[0] == 7);
    static_assert(std::is_const_v<std::remove_reference_t<decltype(bcs)>>,
                  "a const store gives a const sparse set");
  }

  // the walker on the shapes a column may take: only the refs of the moved
  // kind change, whatever the container
  void testRemapWalker() {
    namespace storage = siconos::vkernel::storage;
    namespace config = siconos::vkernel::config;
    using ball_ref = storage::index<config::ball, std::size_t>;
    using inter_ref = storage::index<config::interaction, std::size_t>;

    static_assert(storage::refers_to_v<std::vector<std::variant<ball_ref, inter_ref>>,
                                       config::interaction>);
    static_assert(!storage::holds_refs_v<std::vector<double>>);

    const auto kind = storage::kind_id<config::ball>::get();
    const auto r = storage::ref_remap::swap_and_pop(0, 2);  // 0 dies, 2 -> 0

    std::vector<std::variant<ball_ref, inter_ref>> poly{ball_ref(2), inter_ref(2),
                                                        ball_ref(0)};
    storage::refs_detail::remap<storage::kind_id>(poly, kind, r);
    CPPUNIT_ASSERT_MESSAGE("variant: ball ref 2 -> 0", std::get<0>(poly[0]).value() == 0);
    CPPUNIT_ASSERT_MESSAGE("variant: interaction ref 2 untouched",
                           std::get<1>(poly[1]).value() == 2);
    CPPUNIT_ASSERT_MESSAGE("variant: ball ref 0 dead", !std::get<0>(poly[2]).valid());

    std::map<int, ball_ref> m{{7, ball_ref(2)}, {8, ball_ref(1)}};
    storage::refs_detail::remap<storage::kind_id>(m, kind, r);
    CPPUNIT_ASSERT_MESSAGE("map values", m[7].value() == 0 && m[8].value() == 1);

    std::tuple<ball_ref, std::vector<ball_ref>> t{ball_ref(2), {ball_ref(0), ball_ref(2)}};
    storage::refs_detail::remap<storage::kind_id>(t, kind, r);
    CPPUNIT_ASSERT_MESSAGE("tuple and vector", std::get<0>(t).value() == 0 &&
                                                   !std::get<1>(t)[0].valid() &&
                                                   std::get<1>(t)[1].value() == 0);

    storage::sparse_set<std::size_t, ball_ref> set;
    set.emplace(5, ball_ref(2));
    storage::refs_detail::remap<storage::kind_id>(set, kind, r);
    CPPUNIT_ASSERT_MESSAGE("sparse set values", set.at(5).value() == 0);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(RegistryRemoveTest);

}  // namespace
