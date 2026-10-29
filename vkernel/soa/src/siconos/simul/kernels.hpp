#pragma once

#include "siconos/simul/index_sets.hpp"
#include "siconos/storage/kernels.hpp"

namespace siconos::simul::kernels {

/* The kernels of the simulation (see storage/kernels.hpp): the placeholders
   resolved by the element that runs a kernel, and the iteration spaces of the
   simulation. */

using namespace storage::kernels;

// placeholders: Ctx is the context of the element that runs the kernel
struct system {  // its system type
  template <typename Ctx>
  using resolved_in = typename Ctx::system;
};
struct element {  // the element, declarer of its own properties
  template <typename Ctx>
  using resolved_in = typename Ctx::element;
};
struct topology {  // the topology, declarer of the shared properties
  template <typename Ctx>
  using resolved_in = typename Ctx::topology;
};

// the rank of an involved system among the involved systems
struct rank : space_value {};

// the involved systems of a step, with their rank (see index_sets.hpp: the
// involved set, or the involved and index properties without index sets)
template <typename Item>
struct involved {
  using item = Item;
  template <typename Ctx>
  static void for_each(auto& data, auto step, std::size_t, auto&& f) {
    simul::for_each_involved<resolve_t<Ctx, Item>, resolve_t<Ctx, topology>>(
        data, step, [&f](auto i, auto r) { f(i, r); });
  }
};

}  // namespace siconos::simul::kernels
