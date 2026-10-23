#pragma once

#include <concepts>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/polytope.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* The collision shape of a dynamical system, declared in one of two ways:

   - static: one attached "shape" item_ref per ds type,
       storage::attached<ds, symbol<"shape">, some::item_ref<disk_shape>>
     The shape kind is part of the ds type and every decision is compile time.

   - per kind (ECS sparse sets): one dynamic property per shape kind,
       storage::dynamic_storage<storage::attached<ds, symbol<"disk_shape">,
                                                  some::item_ref<disk_shape>>>
       storage::dynamic_storage<storage::attached<ds, symbol<"polytope_shape">,
                                                  some::item_ref<polytope_shape>>>
     A body has the kind whose set holds it, so a single rigid ds type
     carries disks and polytopes. The sets are registry sparse columns that
     storage::remove keeps consistent.

   The helpers below give the same answers for both declarations. Per-kind
   lookups cost a membership test: callers on a hot path resolve the kind once
   (at pair activation) and keep what they need. */

template <typename Shape>
struct shape_kind;

template <>
struct shape_kind<shape::disk> {
  static constexpr auto name = string_literal("disk_shape");
};

template <>
struct shape_kind<shape::polytope> {
  static constexpr auto name = string_literal("polytope_shape");
};

namespace detail {

template <typename DS>
using ds_item_t = typename std::decay_t<DS>::type;

template <typename DS>
using ds_data_t = std::decay_t<decltype(std::declval<DS&>().data())>;

template <typename DS>
static constexpr bool has_static_shape_v =
    storage::has_property_v<ds_item_t<DS>, ds_data_t<DS>, "shape">;

template <typename Shape, typename DS>
static constexpr bool has_kind_v =
    storage::is_dynamic_storage_v<ds_item_t<DS>, ds_data_t<DS>, shape_kind<Shape>::name,
                                  storage::config_declared>;

template <typename DS>
constexpr auto static_shape_type() {
  using shape_t = typename std::decay_t<decltype(storage::make_handle(
      std::declval<ds_data_t<DS>&>(),
      storage::prop<"shape", storage::config_declared>(std::declval<DS&>())))>::type;
  return std::type_identity<shape_t>{};
}

template <typename Item, typename Shape, typename Data>
decltype(auto) kind_set(Data& data) {
  return storage::prop_dynamic_memory<Item, shape_kind<Shape>::name, storage::config_declared>(
      data);
}

}  // namespace detail

/// DS carries its shape in per-kind sparse sets (at least one kind declared)
template <typename DS>
static constexpr bool per_kind_shapes_v =
    detail::has_kind_v<shape::disk, DS> || detail::has_kind_v<shape::polytope, DS>;

/// every DS item has a Shape (static declaration of a Shape-derived type)
template <typename Shape, typename DS>
constexpr bool always_shape() {
  if constexpr (detail::has_static_shape_v<DS>) {
    return std::derived_from<typename decltype(detail::static_shape_type<DS>())::type, Shape>;
  } else {
    return false;
  }
}

/// some DS item may have a Shape (static match, or a Shape kind is declared)
template <typename Shape, typename DS>
constexpr bool may_have_shape() {
  return always_shape<Shape, DS>() || detail::has_kind_v<Shape, DS>;
}

/// this ds has a Shape. Compile time for a static declaration, a membership
/// test of the Shape set otherwise.
template <typename Shape, typename DS>
bool has_shape(DS& ds) {
  if constexpr (always_shape<Shape, DS>()) {
    return true;
  } else if constexpr (detail::has_kind_v<Shape, DS>) {
    return detail::kind_set<detail::ds_item_t<DS>, Shape>(ds.data()).contains(
        ds.index().value());
  } else {
    return false;
  }
}

/// handle on the Shape of this ds. Throws if the ds has no Shape.
template <typename Shape, typename DS>
auto shape_of(DS& ds) {
  static_assert(may_have_shape<Shape, DS>(), "no ds of this type can have this shape kind");
  auto& data = ds.data();
  if constexpr (always_shape<Shape, DS>()) {
    return storage::make_handle(data, storage::prop<"shape", storage::config_declared>(ds));
  } else {
    auto* ref =
        detail::kind_set<detail::ds_item_t<DS>, Shape>(data).try_get(ds.index().value());
    if (!ref) {
      throw std::runtime_error("collision::shape_of: the body has no shape of this kind");
    }
    return storage::make_handle(data, *ref);
  }
}

/// give this ds the shape sh. With per-kind sets a body has exactly one
/// shape: setting a second kind throws.
template <typename DS, typename SH>
void set_shape(DS& ds, SH sh) {
  using shape_item_t = typename std::decay_t<SH>::type;
  if constexpr (detail::has_static_shape_v<DS>) {
    storage::prop<"shape", storage::config_declared>(ds) = sh;
  } else {
    using kind_t = std::conditional_t<std::derived_from<shape_item_t, shape::polytope>,
                                      shape::polytope, shape::disk>;
    using other_t =
        std::conditional_t<std::same_as<kind_t, shape::disk>, shape::polytope, shape::disk>;
    static_assert(std::derived_from<shape_item_t, kind_t> && detail::has_kind_v<kind_t, DS>,
                  "this ds type declares no set for this shape kind");
    if (has_shape<other_t>(ds)) {
      throw std::runtime_error(
          "collision::set_shape: the body already has a shape of another kind");
    }
    detail::kind_set<detail::ds_item_t<DS>, kind_t>(ds.data()).get_or_create(
        ds.index().value()) = sh;
  }
}

/// the sparse set of the Shape kind of Item: dense_keys_ are the body
/// indices, dense_ the shape refs. For per-kind loops (e.g. every polytope
/// body), without visiting the other bodies.
template <typename Item, typename Shape, typename Data>
decltype(auto) bodies_with(Data& data) {
  return detail::kind_set<Item, Shape>(data);
}

/// radius of the disk shape of ds. Throws when ds has no disk shape, also
/// when no ds of this type can have one: a disk relation is then an
/// alternative of an interaction's relation variant that std::visit
/// instantiates but that is never bound to such a ds. A polytope's radius()
/// is its circumradius, so reading it here would be a silent wrong gap.
template <typename DS>
auto disk_radius(DS& ds) {
  using scalar = typename decltype(ds.env())::scalar;
  if constexpr (may_have_shape<shape::disk, DS>()) {
    return static_cast<scalar>(shape_of<shape::disk>(ds).radius());
  } else {
    throw std::runtime_error(
        "collision::disk_radius: the body has no disk shape; a disk relation was "
        "bound to a body of another shape kind");
    return scalar{};
  }
}

}  // namespace siconos::collision
