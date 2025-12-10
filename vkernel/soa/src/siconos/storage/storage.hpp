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

static auto apply_fun = []<typename Item, typename SomeFun>(
                            auto& data, Item, SomeFun&& some_fun) {
  using item_t = Item;
  using info_t = get_info_t<decltype(data)>;
  using all_keeps_t = decltype(all_properties_as<property::keep>(data));

  using indice = typename info_t::env::indice;

  auto attrs = ground::tuple_unique(
      concat(attributes(item_t{}), attached_storages(item_t{}, data)));

  if constexpr (ground::size(attrs) > ground::size_c<0>) {
    ground::for_each(attrs, [&data, &some_fun]<match::attribute A>(A) {
      return ground::for_each(ground::range<memory_size<A, all_keeps_t>()>,
                              [&data, &some_fun](indice step) {
                                static_cast<SomeFun&&>(some_fun)(memory(
                                    step, ground::get<A>(data.store())));
                              });
    });
  }
};

template <match::item T>
static constexpr void for_each_attribute(T)
{
  return ground::compose(ground::for_each, attributes)(T{});
};

using pattern::attr_t;
using pattern::wrap;

// composite view of attributes and interfaces of 2 items
// properties attached to items are lost
template <typename Item1, typename Item2>
struct composite_item : item<> {
  using items = gather<Item1, Item2>;

  using attributes =
      decltype(concat(std::declval<typename Item1::attributes>(),
                      std::declval<typename Item2::attributes>()));

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    struct composite_holder : Item1::template interface<Handle>,
                              Item2::template interface<Handle> {
      composite_holder(Handle* handlep) : handlep(handlep){};
      // Assignment from Item1 handle (e.g., translated_disk_shape =
      // disk_shape)
      template <typename I, typename R, typename D>
      const composite_holder& operator=(
          const storage::handle<storage::handle_base, I, R, D>& other)
      {
        if constexpr (std::derived_from<Item1, I> ||
                      std::derived_from<Item2, I>) {
          auto& data = handlep->data();

          // Copy all Item1 storages (attributes only)
          Item1 item1{};
          auto storages_to_copy = flatten(all_attributes(item1));

          ground::for_each(
              storages_to_copy,
              [this, &other, &data]<match::attribute Storage>(Storage) {
                ground::get<Storage>(data.store())[(*handlep).get()] =
                  ground::get<Storage>(data.store())[other.index().value()];
              });
          return *this;
        }
        else {
          []<bool flag = false>() {
            static_assert(
                flag, "cannot copy this handle in this composite handle!");
          }();
        }
      }

      Handle* handlep;
    };

    composite_holder composite() { return composite_holder(self()); }
  };
};

}  // namespace siconos::storage
