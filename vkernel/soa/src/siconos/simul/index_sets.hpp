#pragma once

#include <concepts>
#include <type_traits>

#include "siconos/simul/simul_head.hpp"
#include "siconos/storage/sparse_set.hpp"

namespace siconos::simul {

/* The active interactions and the involved systems of a step as sets (ECS
   tags).

   compute_active_interactions decides which interactions are active and
   gives the involved systems their index in the assembled systems. With
   the index sets, it also writes them in one sparse set per item type,
   declared by the topology:

     dynamic_storage<attached<Inter, active_set_tag, indice>>
       interaction index -> its row among the active interactions of Inter
     dynamic_storage<attached<System, involved_set_tag, indice>>
       system index -> its index among the involved systems

   and the assembly loops go over the dense side of these sets
   (for_each_active, for_each_involved) instead of every interaction or
   system. The "activation", "involved" and "index" properties are still
   written: the space filter marks the interactions to keep with
   "activation", io reads them. */

/// tag of the set of the active interactions of a type
struct active_set_tag {};

/// tag of the set of the involved systems of a type
struct involved_set_tag {};

template <typename Item, typename Tag>
using index_set_property = storage::without_binding<
    storage::dynamic_storage<storage::attached<Item, Tag, some::indice>>>;

/// the index set property of Item, none for an absent item type (empty_item)
template <typename Item, typename Tag>
constexpr auto index_set_properties() {
  if constexpr (std::derived_from<Item, empty_item>) {
    return gather<>{};
  } else {
    return gather<index_set_property<Item, Tag>>{};
  }
}

template <typename Item, typename Tag>
using index_set_properties_t = decltype(index_set_properties<Item, Tag>());

/// the topology Decl declares the index sets
template <typename Decl>
static constexpr bool with_index_sets = requires { requires Decl::index_sets; };

/// the set of Item tagged by Tag
template <typename Item, typename Tag, typename Decl>
decltype(auto) index_set(auto& data) {
  return storage::prop_dynamic_memory_tagged<Item, Tag, Decl>(data);
}

/// f(i, row) for every active interaction i of type Inter, row being its rank
/// among the active interactions of Inter
template <typename Inter, typename Decl, typename F>
void for_each_active(auto& data, auto step, F&& f) {
  if constexpr (with_index_sets<Decl>) {
    storage::for_each_entry(index_set<Inter, active_set_tag, Decl>(data), f);
  } else {
    auto& activations = storage::prop_values<Inter, "activation", Decl>(data, step);
    std::size_t row = 0;
    for (std::size_t i = 0; i < std::size(activations); ++i) {
      if (activations[i]) {
        f(i, row++);
      }
    }
  }
}

/// f(i, index) for every involved system i of type System, index being its
/// index among the involved systems
template <typename System, typename Decl, typename F>
void for_each_involved(auto& data, auto step, F&& f) {
  if constexpr (with_index_sets<Decl>) {
    storage::for_each_entry(index_set<System, involved_set_tag, Decl>(data), f);
  } else {
    auto& involveds = storage::prop_values<System, "involved", Decl>(data, step);
    auto& indices = storage::prop_values<System, "index", Decl>(data, step);
    for (std::size_t i = 0; i < std::size(involveds); ++i) {
      if (involveds[i]) {
        f(i, indices[i]);
      }
    }
  }
}

}  // namespace siconos::simul
