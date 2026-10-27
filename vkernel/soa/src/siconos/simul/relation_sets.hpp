#pragma once

#include <cassert>
#include <utility>
#include <variant>

#include "siconos/simul/simul_head.hpp"
#include "siconos/storage/sparse_set.hpp"

namespace siconos::simul {

/* Per relation type sets of interactions (ECS tags).

   An interaction has one relation among the alternatives of its relation
   variant. With the relation sets, its relation type is also given by the
   membership of the interaction in one sparse set per relation type R,
   attached to the interaction type and declared by the topology:

     dynamic_storage<attached<Inter, relation_set_tag<R>, item_ref<R>>>

   with the interaction index as key and the ref of its relation as value.
   The loops of the integrator go over the relation types at compile time and
   over the dense side of each set at run time (for_each_relation): no visit
   per interaction, and what depends on the relation type (one or two ds,
   linear or not) is decided at compile time.

   The sets are kept by the stages that bind the relations: bind_relation
   (topology::bind) is the only way to set a relation, the relation()
   accessor of an interaction is read only. It moves the interaction from
   the set of its previous relation type to the set of the new one (the
   space filter at creation, the cxd narrowphase at its rebind), and
   storage::remove re-keys the sets as it moves the last interaction.
   sync_relation_sets, once per step, checks that the sets hold every
   interaction: an interaction created and not bound (or a relation written
   through storage::attr<"relation">) is a bug, an assertion in debug, and
   the sets are rebuilt from the variants in release. */

/// tag of the set of the interactions whose relation is an R
template <typename R>
struct relation_set_tag {};

template <typename Inter, typename R>
using relation_set_property = storage::without_binding<storage::dynamic_storage<
    storage::attached<Inter, relation_set_tag<R>, some::item_ref<R>>>>;

/// the relation set properties of an interaction type, none for an
/// interaction type without relations (empty_item)
template <typename Inter>
constexpr auto relation_set_properties() {
  if constexpr (requires { typename Inter::relations; }) {
    return []<typename... Rs>(gather<Rs...>) {
      return gather<relation_set_property<Inter, Rs>...>{};
    }(typename Inter::relations{});
  } else {
    return gather<>{};
  }
}

template <typename Inter>
using relation_set_properties_t = decltype(relation_set_properties<Inter>());

/// the set of the interactions of type Inter whose relation is an R
template <typename Inter, typename R, typename Decl>
decltype(auto) relation_set(auto& data) {
  return storage::prop_dynamic_memory_tagged<Inter, relation_set_tag<R>, Decl>(data);
}

/// rebuild the relation sets of Inter from the relation variants
template <typename Inter, typename Decl>
void classify_relations(auto& data, auto step) {
  using relations_t = typename Inter::relations;

  // the sets, looked up once
  auto sets = []<typename... Rs>(auto& data, gather<Rs...>) {
    return std::make_tuple(&relation_set<Inter, Rs, Decl>(data)...);
  }(data, relations_t{});

  std::apply([](auto*... set) { (set->clear(), ...); }, sets);

  auto& relations = storage::attr_values<attr_t<Inter, "relation">>(data, step);
  for (std::size_t i = 0; i < std::size(relations); ++i) {
    std::visit(
        [&sets, i](auto& ref) {
          using R = typename std::decay_t<decltype(ref)>::type;
          using set_t = std::remove_reference_t<decltype(relation_set<Inter, R, Decl>(
              std::declval<std::decay_t<decltype(data)>&>()))>;
          std::get<set_t*>(sets)->emplace(i, ref);
        },
        relations[i]);
  }
}

/// the relation sets of Inter hold every interaction: an assertion in debug,
/// a rebuild from the relation variants in release
template <typename Inter, typename Decl>
void sync_relation_sets(auto& data, auto step) {
  std::size_t n = 0;
  mp::for_each(typename Inter::relations{},
               [&data, &n]<typename R>(R) { n += relation_set<Inter, R, Decl>(data).size(); });
  if (n != std::size(storage::attr_values<attr_t<Inter, "relation">>(data, step))) {
    assert(false && "an interaction is not bound: set its relation with topology::bind");
    classify_relations<Inter, Decl>(data, step);
  }
}

/// the relation of the interaction inter is rel (a handle or a ref on a
/// relation item): its relation variant and, with the relation sets of the
/// topology Decl, its membership
template <typename Decl>
void bind_relation(auto& data, auto inter, auto rel) {
  auto& relvar = storage::attr<"relation">(inter);
  if constexpr (requires { requires Decl::relation_sets; }) {
    using Inter = typename std::decay_t<decltype(inter)>::type;
    const auto i = inter.index().value();
    mp::for_each(typename Inter::relations{},
                 [&data, i]<typename R>(R) { relation_set<Inter, R, Decl>(data).erase(i); });
    relvar = rel;
    std::visit(
        [&data, i](auto& ref) {
          using R = typename std::decay_t<decltype(ref)>::type;
          relation_set<Inter, R, Decl>(data).emplace(i, ref);
        },
        relvar);
  } else {
    relvar = rel;
  }
}

/// f(i, rel) for every interaction i of type Inter, rel being a handle on its
/// relation: relation type by relation type, so that f is instantiated and
/// called with one relation handle type per loop
template <typename Inter, typename Decl, typename F>
void for_each_relation(auto& data, F&& f) {
  mp::for_each(typename Inter::relations{}, [&data, &f]<typename R>(R) {
    auto& set = relation_set<Inter, R, Decl>(data);
    for (std::size_t j = 0; j < set.dense_.size(); ++j) {
      auto rel = storage::make_handle(data, set.dense_[j]);
      f(set.dense_keys_[j], rel);
    }
  });
}

}  // namespace siconos::simul
