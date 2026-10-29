#pragma once

#include "siconos/simul/index_sets.hpp"
#include "siconos/simul/relation_sets.hpp"
#include "siconos/storage/kernels.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul::kernels {

/* The kernels of the simulation (see storage/kernels.hpp): the placeholders
   resolved by the element that runs a kernel, and the iteration spaces of the
   simulation. */

using namespace storage::kernels;

// placeholders: Ctx is the context of the element that runs the kernel
struct systems {  // its system type
  static constexpr std::string_view label = "systems";
  template <typename Ctx>
  using resolved_in = typename Ctx::system;
};
struct interactions {  // its interaction type
  static constexpr std::string_view label = "interactions";
  template <typename Ctx>
  using resolved_in = typename Ctx::interaction;
};
struct by_element {  // the element, declarer of its own properties
  template <typename Ctx>
  using resolved_in = typename Ctx::element;
};
struct by_topology {  // the topology, declarer of the shared properties
  template <typename Ctx>
  using resolved_in = typename Ctx::topology;
};

// the rank of an involved system among the involved systems, of an active
// interaction among the active interactions
struct rank : space_value {};

// the relation of an interaction (a handle)
struct relation : space_value {};

// the involved systems of a step, with their rank (see index_sets.hpp: the
// involved set, or the involved and index properties without index sets)
template <typename Item>
struct involved {
  using item = Item;
  static std::string text() { return "involved " + storage::kernels::detail::label<Item>(); }
  template <typename Ctx>
  static void for_each(auto& data, auto step, std::size_t, auto&& f) {
    simul::for_each_involved<resolve_t<Ctx, Item>, resolve_t<Ctx, by_topology>>(
        data, step, [&f](auto i, auto r) { f(i, r); });
  }
};

// the active interactions of a step, with their rank (see index_sets.hpp)
template <typename Item>
struct active {
  using item = Item;
  static std::string text() { return "active " + storage::kernels::detail::label<Item>(); }
  template <typename Ctx>
  static void for_each(auto& data, auto step, std::size_t, auto&& f) {
    simul::for_each_active<resolve_t<Ctx, Item>, resolve_t<Ctx, by_topology>>(
        data, step, [&f](auto i, auto r) { f(i, r); });
  }
};

// every interaction with its relation: relation type by relation type with
// the relation sets (see relation_sets.hpp), through a visit of the relation
// variant of each interaction otherwise. A kernel on the relations runs once
// per relation type.
template <typename Item>
struct relations {
  using item = Item;
  static std::string text() {
    return "relations of " + storage::kernels::detail::label<Item>();
  }
  template <typename Ctx>
  static void for_each(auto& data, auto step, std::size_t, auto&& f) {
    using inter_t = resolve_t<Ctx, Item>;
    using topology_t = resolve_t<Ctx, by_topology>;
    if constexpr (requires { requires topology_t::relation_sets; }) {
      simul::for_each_relation<inter_t, topology_t>(data, f);
    } else {
      auto& relations = storage::attr_values<inter_t, "relation">(data, step);
      for (std::size_t i = 0; i < std::size(relations); ++i) {
        siconos::variant::visit(data, relations[i], [&f, i](auto& rel) { f(i, rel); });
      }
    }
  }
};

}  // namespace siconos::simul::kernels
