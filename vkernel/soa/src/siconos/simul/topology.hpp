#pragma once

#include "siconos/algebra/sparse_factor.hpp"
#include "siconos/simul/index_sets.hpp"
#include "siconos/simul/relation_sets.hpp"
#include "siconos/simul/simul_head.hpp"
#include "siconos/utils/range.hpp"

// default of the RelationSets switch of the topology: the relation type of an
// interaction given by per relation type sets (1) or by its relation variant
// alone (0)
#ifndef VKERNEL_RELATION_SETS
#define VKERNEL_RELATION_SETS 1
#endif

// default of the IndexSets switch of the topology: the active interactions
// and the involved systems of a step also given as sets (1) or by their
// boolean properties alone (0)
#ifndef VKERNEL_INDEX_SETS
#define VKERNEL_INDEX_SETS 1
#endif

namespace siconos::simul {

template <typename FixedDofDynamicalSystem, typename FixedDofInteraction,
          typename DynamicDofDynamicalSystem = empty_item,
          typename DynamicDofFixedDofInteraction = empty_item,
          typename DynamicDofDynamicDofInteraction = empty_item,
          bool BcVelocitiesDynamic = true, bool RelationSets = VKERNEL_RELATION_SETS,
          bool IndexSets = VKERNEL_INDEX_SETS>
struct topology : item {
  // fixed dof case
  using dof = some::indice_parameter<"dof">;

  using fixed_dof_system = FixedDofDynamicalSystem;
  using dynamic_dof_system = DynamicDofDynamicalSystem;
  using fixed_dof_interaction = FixedDofInteraction;
  using dynamic_dof_fixed_dof_interaction = DynamicDofFixedDofInteraction;
  using dynamic_dof_dynamic_dof_interaction = DynamicDofDynamicDofInteraction;

  using items = gather<fixed_dof_system, dynamic_dof_system, dynamic_dof_fixed_dof_interaction,
                       dynamic_dof_dynamic_dof_interaction>;

  using nslaw = typename fixed_dof_interaction::nslaw;
  using nslaw_size = some::indice_value<nslaw::size>;

  using fsystem = fixed_dof_system;
  using dsystem = dynamic_dof_system;

  using finteraction = fixed_dof_interaction;
  using dfinteraction = dynamic_dof_fixed_dof_interaction;
  using ddinteraction = dynamic_dof_dynamic_dof_interaction;

  struct attributes {
    some::map<some::indice, some::item_ref<fsystem>> system_id;
  };

  using bc_velocities_0_fsys_attr = storage::attached<fsystem, symbol<"bc_velocities_0">,
                                                      some::unbounded_vector<some::indice>>;
  using bc_velocities_0_dsys_attr =
      storage::attached<dsystem, symbol<"bc_velocities_0">,
                        some::unbounded_collection<some::indice>>;

  using bc_velocities_0_fsys =
      std::conditional_t<BcVelocitiesDynamic,
                         storage::dynamic_storage<bc_velocities_0_fsys_attr>,
                         bc_velocities_0_fsys_attr>;
  using bc_velocities_0_dsys =
      std::conditional_t<BcVelocitiesDynamic,
                         storage::dynamic_storage<bc_velocities_0_dsys_attr>,
                         bc_velocities_0_dsys_attr>;

  using base_properties =
      gather<storage::attached<fsystem, symbol<"involved">, some::boolean>,
             storage::attached<fsystem, symbol<"index">, some::indice>,
             storage::attached<fsystem, symbol<"id">, some::indice>,
             storage::attached<
                 fsystem, symbol<"p0">,
                 some::array<some::vector<some::scalar, dof>, std::integral_constant<int, 2>>>,

             bc_velocities_0_fsys,

             storage::attached<dsystem, symbol<"q0">, some::unbounded_vector<some::scalar>>,
             storage::attached<dsystem, symbol<"involved">, some::boolean>,
             storage::attached<dsystem, symbol<"index">, some::indice>,
             storage::attached<dsystem, symbol<"id">, some::indice>,
             storage::without_binding<storage::attached<
                 dsystem, symbol<"p0">,
                 some::array<some::unbounded_vector<
                                 some::vector<some::scalar, std::integral_constant<int, 1>>>,
                             std::integral_constant<int, 2>>>>,

             bc_velocities_0_dsys,

             // factorizations of M + h^2 theta^2 K and M, built once (see
             // algebra::fem_operator)
             storage::without_binding<
                 storage::attached<dsystem, symbol<"fem_operator">,
                                   some::specific<std::shared_ptr<algebra::fem_operator>>>>,

             storage::attached<finteraction, symbol<"nds">, some::indice>,
             storage::attached<finteraction, symbol<"ds1">, some::item_ref<fsystem>>,
             storage::attached<finteraction, symbol<"ds2">, some::item_ref<fsystem>>,
             storage::attached<finteraction, symbol<"activation">, some::boolean>,

             storage::attached<dfinteraction, symbol<"nds">, some::indice>,
             storage::attached<dfinteraction, symbol<"ds1">, some::item_ref<fsystem>>,
             storage::attached<dfinteraction, symbol<"ds2">, some::item_ref<dsystem>>,
             storage::attached<dfinteraction, symbol<"activation">, some::boolean>,

             storage::attached<ddinteraction, symbol<"nds">, some::indice>,
             storage::attached<ddinteraction, symbol<"ds1">, some::item_ref<dsystem>>,
             storage::attached<ddinteraction, symbol<"ds2">, some::item_ref<dsystem>>,
             storage::attached<ddinteraction, symbol<"activation">, some::boolean>>;

  // the relation type of an interaction as the membership of one set per
  // relation type (see relation_sets.hpp), or the relation variant alone
  static constexpr bool relation_sets = RelationSets;

  using relation_sets_properties =
      std::conditional_t<RelationSets,
                         decltype(concat(
                             concat(simul::relation_set_properties_t<finteraction>{},
                                    simul::relation_set_properties_t<dfinteraction>{}),
                             simul::relation_set_properties_t<ddinteraction>{})),
                         gather<>>;

  // the active interactions and the involved systems of a step as sets (see
  // index_sets.hpp), or their boolean properties alone
  static constexpr bool index_sets = IndexSets;

  using index_sets_properties = std::conditional_t<
      IndexSets,
      decltype(concat(
          concat(concat(simul::index_set_properties_t<finteraction, active_set_tag>{},
                        simul::index_set_properties_t<dfinteraction, active_set_tag>{}),
                 simul::index_set_properties_t<ddinteraction, active_set_tag>{}),
          concat(simul::index_set_properties_t<fsystem, involved_set_tag>{},
                 simul::index_set_properties_t<dsystem, involved_set_tag>{}))),
      gather<>>;

  using properties = decltype(concat(concat(base_properties{}, relation_sets_properties{}),
                                     index_sets_properties{}));

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    /* self rigid */
    template <match::handle<fsystem> Hds>
    decltype(auto) link(Hds ds) {
      auto& data = self()->data();

      auto inter = storage::add<finteraction>(data);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      prop<"nds", topology>(inter) = 1;
      prop<"ds1", topology>(inter) = ds;
      prop<"ds2", topology>(inter) = ds;

      return inter;
    };

    /* rigid <-> rigid */
    template <match::handle<fsystem> Hds>
    decltype(auto) link(Hds ds1, Hds ds2) {
      auto& data = self()->data();

      auto inter = storage::add<finteraction>(data);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      prop<"nds", topology>(inter) = 2;
      prop<"ds1", topology>(inter) = ds1;
      prop<"ds2", topology>(inter) = ds2;

      return inter;
    };

    /* rigid <-> fem */
    template <match::handle<fsystem> Hfds, match::handle<dsystem> Hdds>
    decltype(auto) link(Hfds ds1, Hdds ds2) {
      auto& data = self()->data();

      auto inter = storage::add<dfinteraction>(data);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      prop<"nds", topology>(inter) = 2;

      prop<"ds1", topology>(inter) = ds1;
      prop<"ds2", topology>(inter) = ds2;

      return inter;
    };

    /* self fem */
    template <match::handle<dsystem> Hdds>
    decltype(auto) link(Hdds ds) {
      auto& data = self()->data();

      auto inter = storage::add<ddinteraction>(data);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      prop<"nds", topology>(inter) = 1;
      prop<"ds1", topology>(inter) = ds;
      prop<"ds2", topology>(inter) = ds;

      return inter;
    };

    /* the relation of the interaction inter is rel (a handle or a ref on a
       relation item): the only way to set it, as it also keeps the relation
       sets (see relation_sets.hpp) */
    void bind(auto inter, auto rel) {
      simul::bind_relation<topology>(self()->data(), inter, rel);
    }

    /* reorder the interactions of each type by relation type (see
       simul::sort_interactions), so that the loops over a relation type go
       through a contiguous range of the interaction storages. To be called
       between two steps, e.g. every few steps as the space filter adds and
       removes interactions. It changes the indices of the interactions,
       hence the order of the rows of the nonsmooth problem (and the
       iterates of an order dependent solver such as NSGS): no handle on an
       interaction may be held across the call. */
    void sort_interactions(auto step) {
      auto& data = self()->data();
      mp::for_each(gather<finteraction, dfinteraction, ddinteraction>{},
                   [&data, step]<typename Inter>(Inter) {
                     if constexpr (!std::derived_from<Inter, empty_item>) {
                       simul::sort_interactions<Inter, topology>(data, step);
                     }
                   });
    }

    void set_dynamical_system_id(auto hds, auto id) {
      attr<"system_id">(*self())[id] = hds.index();
    }

    decltype(auto) dynamical_system(auto id) {
      return storage::make_handle(self()->data(), attr<"system_id">(*self())[id]);
    }

    auto methods() {
      using data_t = std::decay_t<decltype(self()->data())>;
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using hds_t = storage::handle<storage::handle_base, fsystem, indice, data_t>;

      return collect(
          method("set_dynamical_system_id",
                 &interface<Handle>::set_dynamical_system_id<hds_t, indice>),
          method("dynamical_system", &interface<Handle>::dynamical_system<indice>),
          method("sort_interactions", &interface<Handle>::sort_interactions<indice>));
    }
  };
};
}  // namespace siconos::simul
