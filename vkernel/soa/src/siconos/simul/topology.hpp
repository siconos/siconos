#pragma once

#include "siconos/simul/simul_head.hpp"
#include "siconos/utils/print.hpp"
#include "siconos/utils/range.hpp"

namespace siconos::simul {

template <typename FixedDofDynamicalSystem, typename FixedDofInteraction,
          typename DynamicDofDynamicalSystem = empty_item,
          typename DynamicDofFixedDofInteraction = empty_item,
          typename DynamicDofDynamicDofInteraction = empty_item>
struct topology : item {
  // fixed dof case
  using dof = some::indice_parameter<"dof">;

  using fixed_dof_system = FixedDofDynamicalSystem;
  using dynamic_dof_system = DynamicDofDynamicalSystem;
  using fixed_dof_interaction = FixedDofInteraction;
  using dynamic_dof_fixed_dof_interaction = DynamicDofFixedDofInteraction;
  using dynamic_dof_dynamic_dof_interaction = DynamicDofDynamicDofInteraction;

  using items = gather<fixed_dof_system, dynamic_dof_system,
                       dynamic_dof_fixed_dof_interaction,
                       dynamic_dof_dynamic_dof_interaction>;

  using nslaw = typename fixed_dof_interaction::nslaw;
  using nslaw_size = some::indice_value<nslaw::size>;

  using fsystem = fixed_dof_system;
  using dsystem = dynamic_dof_system;

  using finteraction = fixed_dof_interaction;
  using dfinteraction = dynamic_dof_fixed_dof_interaction;
  using ddinteraction = dynamic_dof_dynamic_dof_interaction;

  using attributes =
      gather<attribute<"system_id",
                       some::map<some::indice, some::item_ref<fsystem>>>>;

  using properties = gather<
      storage::attached<fsystem, symbol<"involved">, some::boolean>,
      storage::attached<fsystem, symbol<"index">, some::indice>,
      storage::attached<fsystem, symbol<"id">, some::indice>,
      storage::attached<fsystem, symbol<"p0">,
                        some::array<some::vector<some::scalar, dof>,
                                    std::integral_constant<int, 2>>>,

      storage::attached<fsystem, symbol<"bc_velocities_0">,
                        some::unbounded_collection<some::indice>>,

      storage::attached<dsystem, symbol<"involved">, some::boolean>,
      storage::attached<dsystem, symbol<"index">, some::indice>,
      storage::attached<dsystem, symbol<"id">, some::indice>,
      storage::attached<
          dsystem, symbol<"p0">,
          some::array<some::unbounded_vector<some::vector<
                          some::scalar, std::integral_constant<int, 1>>>,
                      std::integral_constant<int, 2>>>,

      storage::attached<dsystem, symbol<"bc_velocities_0">,
                        some::unbounded_collection<some::indice>>,

      storage::attached<finteraction, symbol<"nds">, some::indice>,
      storage::attached<finteraction, symbol<"ds1">, some::item_ref<fsystem>>,
      storage::attached<finteraction, symbol<"ds2">, some::item_ref<fsystem>>,
      storage::attached<finteraction, symbol<"ydot_backup">,
                        some::vector<some::scalar, nslaw_size>>,
      storage::attached<finteraction, symbol<"activation">, some::boolean>,

      storage::attached<dfinteraction, symbol<"nds">, some::indice>,
      storage::attached<dfinteraction, symbol<"ds1">,
                        some::item_ref<fsystem>>,
      storage::attached<dfinteraction, symbol<"ds2">,
                        some::item_ref<dsystem>>,
      storage::attached<dfinteraction, symbol<"ydot_backup">,
                        some::vector<some::scalar, nslaw_size>>,
      storage::attached<dfinteraction, symbol<"activation">, some::boolean>,

      storage::attached<ddinteraction, symbol<"nds">, some::indice>,
      storage::attached<ddinteraction, symbol<"ds1">,
                        some::item_ref<dsystem>>,
      storage::attached<ddinteraction, symbol<"ds2">,
                        some::item_ref<dsystem>>,
      storage::attached<ddinteraction, symbol<"ydot_backup">,
                        some::vector<some::scalar, nslaw_size>>,
      storage::attached<ddinteraction, symbol<"activation">, some::boolean>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    /* self rigid */
    template <match::handle<fsystem> Hds>
    decltype(auto) link(Hds ds)
    {
      auto &data = self()->data();

      auto inter = storage::add<finteraction>(data);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      prop<"nds">(inter) = 1;
      prop<"ds1">(inter) = ds;
      prop<"ds2">(inter) = ds;

      return inter;
    };

    /* rigid <-> rigid */
    template <match::handle<fsystem> Hds>
    decltype(auto) link(Hds ds1, Hds ds2)
    {
      auto &data = self()->data();

      auto inter = storage::add<finteraction>(data);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      prop<"nds">(inter) = 2;
      prop<"ds1">(inter) = ds1;
      prop<"ds2">(inter) = ds2;

      return inter;
    };

    /* rigid <-> fem */
    template <match::handle<fsystem> Hfds, match::handle<dsystem> Hdds>
    decltype(auto) link(Hfds ds1, Hdds ds2)
    {
      auto &data = self()->data();

      auto inter = storage::add<dfinteraction>(data);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      prop<"nds">(inter) = 2;

      prop<"ds1">(inter) = ds1;
      prop<"ds2">(inter) = ds2;

      return inter;
    };

    /* self fem */
    template <match::handle<dsystem> Hdds>
    decltype(auto) link(Hdds ds)
    {
      auto &data = self()->data();

      auto inter = storage::add<ddinteraction>(data);

      attr<"y">(inter).setZero();
      attr<"ydot">(inter).setZero();
      attr<"lambda">(inter).setZero();

      prop<"nds">(inter) = 1;
      prop<"ds1">(inter) = ds;
      prop<"ds2">(inter) = ds;

      return inter;
    };

    void set_dynamical_system_id(auto hds, auto id)
    {
      attr<"system_id">(*self())[id] = hds.index();
    }

    decltype(auto) dynamical_system(auto id)
    {
      return storage::make_handle(self()->data(),
                                  attr<"system_id">(*self())[id]);
    }

    auto methods()
    {
      using data_t = std::decay_t<decltype(self()->data())>;
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using hds_t =
          storage::handle<storage::handle_base, fsystem, indice, data_t>;

      return collect(
          method("set_dynamical_system_id",
                 &interface<Handle>::set_dynamical_system_id<hds_t, indice>),
          method("dynamical_system",
                 &interface<Handle>::dynamical_system<indice>));
    }
  };
};
}  // namespace siconos::simul
