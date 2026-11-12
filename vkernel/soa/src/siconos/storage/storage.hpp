#pragma once

#include "siconos/storage/add.hpp"
#include "siconos/storage/data_holder.hpp"
#include "siconos/storage/get.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/make.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/remove.hpp"

namespace siconos::storage {

// template <match::item T>
// static constexpr void for_each_attribute(T)
// {
//   return mp::compose(mp::for_each, attributes)(T{});
// };

using pattern::attr_t;
using pattern::wrap;

// composite view of attributes and interfaces of 2 items
// properties attached to items are lost
template <typename Item1, typename Item2>
struct composite_item : item {
  using items = gather<Item1, Item2>;

  using attributes = decltype(concat(std::declval<typename Item1::attributes>(),
                                     std::declval<typename Item2::attributes>()));

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    struct composite_holder : Item1::template interface<Handle>,
                              Item2::template interface<Handle> {
      composite_holder(Handle* handlep) : handlep(handlep) {};
      // Assignment from Item1 handle (e.g., translated_disk_shape =
      // disk_shape)
      template <typename I, typename R, typename D>
      const composite_holder& operator=(
          const storage::handle<storage::handle_base, I, R, D>& other) {
        if constexpr (std::derived_from<Item1, I> || std::derived_from<Item2, I>) {
          auto& data = handlep->data();

          // Copy all Item1 storages (attributes only)
          Item1 item1{};
          auto storages_to_copy = flatten(all_attributes(item1));

          mp::for_each(storages_to_copy,
                       [this, &other, &data]<match::attribute Storage>(Storage) {
                         mp::get<Storage>(data.store())[(*handlep).get()] =
                             mp::get<Storage>(data.store())[other.index().value()];
                       });
          return *this;
        } else {
          []<bool flag = false>() {
            static_assert(flag, "cannot copy this handle in this composite handle!");
          }();
        }
      }

      Handle* handlep;
    };

    composite_holder composite() { return composite_holder(self()); }
  };
};

}  // namespace siconos::storage
