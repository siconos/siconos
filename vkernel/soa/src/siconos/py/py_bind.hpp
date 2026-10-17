#pragma once

// Python bindings of the disks world, one item at a time, so that they can be
// spread over several translation units (see py_bind_part.hpp and main.cpp).

#include <format>
#include <utility>

#include "py_head.hpp"

namespace siconos::python::disks {

using namespace boost::hana::literals;
namespace mp = siconos::storage::mp;

inline data_t make_storage() { return data_t(); };

using disks_info_t = siconos::storage::get_info_t<idata_t>;

using indice_t = typename disks_info_t::template env<config::disk>::indice;

using disks_properties_t = typename disks_info_t::all_properties_t;

using disks_items_t =
    decltype(mp::transform(typename disks_info_t::all_items_t{}, []<match::item I>(I) {
      if constexpr (match::wrap<I>) {
        return typename I::type{};
      } else {
        return I{};
      }
    }));

// the storage::bind<> properties, gathered once
using disks_binds_t =
    decltype(mp::filter(disks_properties_t{}, mp::derive_from<storage::property::bind>));

// T has a bind<> name: what has_property_from<T, property::bind, properties>
// tests, without filtering the whole property list again for every item
template <typename T, typename... Binds>
constexpr bool bound_by(boost::hana::tuple<Binds...>) {
  return (std::derived_from<T, typename Binds::item> || ...);
}

// the items with a storage::bind<> name, in a fixed order shared by all parts
using named_disks_items_t = decltype(mp::tuple_unique(mp::filter(
    disks_items_t{},
    mp::is_a_model<[]<typename T>() { return bound_by<T>(disks_binds_t{}); }>)));

/// py::class_ for the handle (and index, and batch wrapper) of item I, with
/// its attached storages, attributes, methods and add_ functions.
template <match::item I>
void bind_item(py::module_& disks) {
  using indice_i = typename disks_info_t::template env<I>::indice;
  using handle_t = storage::handle<storage::handle_base, I, indice_i, idata_t>;
  using item_t = typename handle_t::type;
  using index_t = typename handle_t::index_t;

  auto pyhandle = [&disks]() {
    auto base_index = py::class_<index_t>(
        disks,
        std::format("index_{}", storage::bind_name<item_t, disks_properties_t>()).c_str());
    if constexpr (match::batch_capable<item_t> && !match::without_batch_bindings<item_t>) {
      return mp::make_tuple(
          base_index,
          py::class_<handle_t>(disks, storage::bind_name<item_t, disks_properties_t>(),
                               base_index),
          mp::type_c<handle_t>,
          py::class_<handles_wrap<std::vector<handle_t>>>(
              disks,
              std::format("multiple_{}", storage::bind_name<item_t, disks_properties_t>())
                  .c_str()));
    }

    else {
      return mp::make_tuple(
          base_index,
          py::class_<handle_t>(disks, storage::bind_name<item_t, disks_properties_t>(),
                               base_index),
          mp::type_c<handle_t>);
    }
  }();

  // attached storage
  if constexpr (!match::without_attached_storages_bindings<item_t>) {
    mp::fold_left(
        decltype(storage::attached_storages(item_t{},
                                            std::declval<idata_t&>())){},  // all attached
                                                                           // storages
        std::ref(pyhandle[1_c]),                                           // initial state
        []<match::attached_storage<item_t> S>(py::class_<handle_t> dc, S s) {
          if constexpr (!match::without_binding<S>) {
            constexpr auto astor_name = storage::attached_storage_name(s);
            using target_type = std::decay_t<decltype(out_formatter(
                handle_t{}, storage::prop<astor_name.str>(handle_t{})))>;

            return dc
                .def(
                    std::format("{}", astor_name.str.value).c_str(),
                    [&astor_name](handle_t& h) -> target_type {
                      return out_formatter(h, storage::prop<astor_name.str>(h));
                    },
                    py::return_value_policy::reference)
                .def(std::format("set_{}", astor_name.str.value).c_str(),
                     [&astor_name](handle_t& h, target_type val) {
                       in_formatter(h, storage::prop<astor_name.str>(h)) =
                           in_formatter(h, val);
                     });
          } else {
            return dc;
          }
        });
  }

  if constexpr (match::batch_capable<item_t> && !match::without_batch_bindings<item_t>) {
    mp::fold_left(
        decltype(storage::attached_storages(item_t{},
                                            std::declval<idata_t&>())){},  // all attached
                                                                           // storages
        std::ref(pyhandle[3_c]),                                           // initial state
        []<match::attached_storage<item_t> S>(
            py::class_<handles_wrap<std::vector<handle_t>>> hw, S s) {
          constexpr auto astor_name = storage::attached_storage_name(s);
          using target_type = std::decay_t<decltype(out_formatter(
              handle_t{}, storage::prop<astor_name.str>(handle_t{})))>;

          return hw
              .def(
                  std::format("multiple_{}", astor_name.str.value).c_str(),
                  [&astor_name](handles_wrap<std::vector<handle_t>>& h_multiple) {
                    return std::ranges::to<std::vector>(
                        h_multiple.handles | view::transform([&astor_name](auto&& h) {
                          return out_formatter(h, storage::prop<astor_name.str>(h));
                        }));
                  },
                  py::return_value_policy::reference)
              .def(std::format("set_{}", astor_name.str.value).c_str(),
                   [&astor_name](handles_wrap<std::vector<handle_t>>& h_multiple,
                                 target_type val) {
                     std::ranges::for_each(
                         h_multiple.handles, [&astor_name, &val](auto&& h) {
                           in_formatter(h, storage::prop<astor_name.str>(h)) =
                               in_formatter(h, val);
                         });
                   })
              .def(std::format("multiple_set_{}", astor_name.str.value).c_str(),
                   [&astor_name](handles_wrap<std::vector<handle_t>>& h_multiple,
                                 std::vector<target_type> vals) {
                     std::ranges::for_each(
                         view::zip(vals, h_multiple.handles), [&astor_name](auto&& pair) {
                           auto& val = std::get<0>(pair);
                           auto& h = std::get<1>(pair);
                           in_formatter(h, storage::prop<astor_name.str>(h)) =
                               in_formatter(h, val);
                         });
                   });
        });
  }

  if constexpr (!match::without_attributes_bindings<typename handle_t::type>) {
    mp::fold_left(
        pattern::attributes(typename handle_t::type{}),  // all attributes
        std::ref(pyhandle[1_c]),                         // initial state
        []<match::attribute A>(py::class_<handle_t> dc, A a) {
          using attr_value_t = std::decay_t<decltype(out_formatter(
              handle_t{}, storage::get<A>(std::declval<idata_t&>(), 0, handle_t{})))>;
          return dc
              .def(
                  std::format("{}", pattern::attribute_name(a)).c_str(),
                  [](handle_t& h) -> attr_value_t {
                    return out_formatter(h, storage::get<A>(h.data(), h));
                  },
                  py::return_value_policy::reference)
              .def(
                  std::format("{}_at_step", pattern::attribute_name(a)).c_str(),
                  [](handle_t& h, indice_t step) -> attr_value_t {
                    return out_formatter(h, storage::get<A>(h.data(), step, h));
                  },
                  py::return_value_policy::reference)

              .def(std::format("set_{}", pattern::attribute_name(a)).c_str(),
                   [](handle_t& h, attr_value_t v) {
                     in_formatter(h, storage::get<A>(h.data(), h)) = in_formatter(h, v);
                   })
              .def(std::format("set_{}_at_step", pattern::attribute_name(a)).c_str(),
                   [](handle_t& h, attr_value_t v, indice_t step) {
                     in_formatter(h, storage::get<A>(h.data(), step, h)) = in_formatter(h, v);
                   });
        });
  }

  if constexpr (match::batch_capable<item_t> && !match::without_batch_bindings<item_t>) {
    mp::fold_left(
        pattern::attributes(typename handle_t::type{}),  // all attributes
        std::ref(pyhandle[3_c]),                         // initial state
        []<match::attribute A>(py::class_<handles_wrap<std::vector<handle_t>>> dc, A a) {
          using attr_value_t = std::decay_t<decltype(out_formatter(
              handle_t{}, storage::get<A>(std::declval<idata_t&>(), 0, handle_t{})))>;
          return dc
              .def("get",
                   [](handles_wrap<std::vector<handle_t>>& h_multiple) {
                     return h_multiple.handles;
                   })
              .def(
                  std::format("multipe_{}", pattern::attribute_name(a)).c_str(),
                  [](handles_wrap<std::vector<handle_t>>& h_multiple) {
                    return std::ranges::to<std::vector>(
                        h_multiple.handles | view::transform([](auto&& h) {
                          return out_formatter(h, storage::get<A>(h.data(), h));
                        }));
                  },
                  py::return_value_policy::reference)
              .def(std::format("set_{}", pattern::attribute_name(a)).c_str(),
                   [](handles_wrap<std::vector<handle_t>>& h_multiple, attr_value_t val) {
                     std::ranges::for_each(h_multiple.handles, [&val](auto&& h) {
                       in_formatter(h, storage::get<A>(h.data(), h)) = in_formatter(h, val);
                     });
                   })
              .def(std::format("multiple_set_{}", pattern::attribute_name(a)).c_str(),
                   [](handles_wrap<std::vector<handle_t>>& h_multiple,
                      std::vector<attr_value_t> vs) {
                     std::ranges::for_each(view::zip(vs, h_multiple.handles),
                                           [](auto&& pair) {
                                             auto& v = std::get<0>(pair);
                                             auto& h = std::get<1>(pair);
                                             in_formatter(h, storage::get<A>(h.data(), h)) =
                                                 in_formatter(h, v);
                                           });
                   });
        });
  }

  mp::fold_left(storage::methods(pyhandle[2_c]),  // all methods
                std::ref(pyhandle[1_c]),          // initial state
                []<typename M>(py::class_<handle_t> dc, M m) {
                  return dc.def(pattern::method_name(m), pattern::method_def(m),
                                py::return_value_policy::reference);
                });

  auto item_name = storage::bind_name<item_t, disks_properties_t>();

  disks.def(
      std::format("add_{}", item_name).c_str(),
      [](data_t& data) { return siconos::storage::add<item_t>(data()); },
      py::return_value_policy::reference, R"pbdoc(
        Add
    )pbdoc");

  if constexpr (match::batch_capable<item_t> && !match::without_batch_bindings<item_t>) {
    disks.def(
        std::format("multiple_add_{}", item_name).c_str(),
        [](data_t& data, indice_t count) {
          auto r = siconos::storage::add<item_t>(data(), count);
          return handles_wrap<std::vector<handle_t>>{std::ranges::to<std::vector<handle_t>>(r)};
        },
        R"pbdoc(
        Add
    )pbdoc");
  }
}

/// Bind the Idx-th named item if it falls in part Part out of Parts.
template <std::size_t Idx, std::size_t Part, std::size_t Parts>
void bind_item_in_part(py::module_& disks) {
  if constexpr (Idx % Parts == Part) {
    bind_item<std::decay_t<decltype(boost::hana::at_c<Idx>(named_disks_items_t{}))>>(disks);
  }
}

/// Bind the named items of part Part: items are dealt round-robin over Parts.
template <std::size_t Part, std::size_t Parts>
void bind_items(py::module_& disks) {
  constexpr std::size_t n = decltype(mp::size(named_disks_items_t{}))::value;
  [&disks]<std::size_t... Idx>(std::index_sequence<Idx...>) {
    (bind_item_in_part<Idx, Part, Parts>(disks), ...);
  }(std::make_index_sequence<n>{});
}

}  // namespace siconos::python::disks
