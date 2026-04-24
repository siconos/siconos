#include "py_head.hpp"

namespace siconos::python::disks {
data_t make_storage() { return data_t(); };
}  // namespace siconos::python::disks

using namespace boost::hana::literals;
namespace mp = siconos::storage::mp;

PYBIND11_MODULE(_nonos, m)
{
  // sub module disks
  auto disks = m.def_submodule("disks");

  disks.doc() = R"pbdoc(
        Nonos m
        -----------

        .. currentmodule:: nonos

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

  disks.def("make_storage", &siconos::python::disks::make_storage,
            py::return_value_policy::reference,
            R"pbdoc(
        Create a new data object for 2D disks simulation
    )pbdoc");

  auto data_class =
      py::class_<siconos::python::disks::data_t>(disks, "data_t");

  using disks_info_t =
      siconos::storage::get_info_t<siconos::python::disks::idata_t>;

  using indice_t =
      typename disks_info_t::template env<config::disks::disk>::indice;

  using disks_properties_t = typename disks_info_t::all_properties_t;

  using disks_items_t = decltype(mp::transform(
      typename disks_info_t::all_items_t{}, []<match::item I>(I) {
        if constexpr (match::wrap<I>) {
          return typename I::type{};
        }
        else {
          return I{};
        }
      }));

  // mp::type_trace<disks_items_t>();
  auto named_disks_items = mp::tuple_unique(mp::filter(
      disks_items_t{}, mp::is_a_model<[]<typename T>() {
        return storage::has_property_from<T, storage::property::bind,
                                          disks_properties_t>();
      }>));

  // mp::type_trace<std::decay_t<decltype(named_disks_items)>>();

  auto disks_handles =
      mp::transform(named_disks_items, []<match::item I>(I item) constexpr {
        using indice_t = typename storage::get_info_t<
            siconos::python::disks::idata_t>::template env<I>::indice;
        using handle_t = storage::handle<storage::handle_base, I, indice_t,
                                         siconos::python::disks::idata_t>;
        return mp::type_c<handle_t>;
      });

  // add corresponding py::class_
  auto pyhandles = mp::transform(disks_handles, [&disks]<typename H>(
                                                    H handle_type_c) {
    using handle_t =
        typename decltype(handle_type_c)::type;  // Extract actual type
    using item_t = typename handle_t::type;
    using index_t = typename handle_t::index_t;
    auto base_index = py::class_<index_t>(
        disks, fmt::format("index_{}",
                           storage::bind_name<item_t, disks_properties_t>())
                   .c_str());
    if constexpr (match::batch_capable<item_t>) {
      return mp::make_tuple(
          base_index,
          py::class_<handle_t>(
              disks, storage::bind_name<item_t, disks_properties_t>(),
              base_index),
          mp::type_c<handle_t>,
          py::class_<handles_wrap<std::vector<handle_t>>>(
              disks,
              fmt::format("multiple_{}",
                          storage::bind_name<item_t, disks_properties_t>())
                  .c_str()));
    }

    else {
      return mp::make_tuple(
          base_index,
          py::class_<handle_t>(
              disks, storage::bind_name<item_t, disks_properties_t>(),
              base_index),
          mp::type_c<handle_t>);
    }
  });

  // attached storage
  mp::for_each(pyhandles, [&disks](auto pyhandle) {
    using handle_t = typename decltype(+pyhandle[2_c])::type;
    using item_t = typename handle_t::type;

    if constexpr (!match::without_attached_storages_bindings<item_t>) {
      mp::fold_left(
          decltype(storage::attached_storages(
              item_t{},
              std::declval<
                  siconos::python::disks::idata_t&>())){},  // all attached
                                                            // storages
          std::ref(pyhandle[1_c]),                          // initial state
          []<match::attached_storage<item_t> S>(py::class_<handle_t> dc,
                                                S s) {
            constexpr auto astor_name = storage::attached_storage_name(s);
            using target_type = std::decay_t<decltype(out_formatter(
                handle_t{}, storage::prop<astor_name.str>(handle_t{})))>;

            return dc
                .def(
                    fmt::format("{}", astor_name.str.value).c_str(),
                    [&astor_name](handle_t& h) -> target_type {
                      return out_formatter(h,
                                           storage::prop<astor_name.str>(h));
                    },
                    py::return_value_policy::reference)
                .def(fmt::format("set_{}", astor_name.str.value).c_str(),
                     [&astor_name](handle_t& h, target_type val) {
                       in_formatter(h, storage::prop<astor_name.str>(h)) =
                           in_formatter(h, val);
                     });
          });
    }

    if constexpr (match::batch_capable<item_t>) {
      mp::fold_left(
          decltype(storage::attached_storages(
              item_t{},
              std::declval<
                  siconos::python::disks::idata_t&>())){},  // all attached
                                                            // storages
          std::ref(pyhandle[3_c]),                          // initial state
          []<match::attached_storage<item_t> S>(
              py::class_<handles_wrap<std::vector<handle_t>>> hw, S s) {
            constexpr auto astor_name = storage::attached_storage_name(s);
            using target_type = std::decay_t<decltype(out_formatter(
                handle_t{}, storage::prop<astor_name.str>(handle_t{})))>;

            return hw
                .def(
                    fmt::format("multiple_{}", astor_name.str.value).c_str(),
                    [&astor_name](
                        handles_wrap<std::vector<handle_t>>& h_multiple) {
                      return ranges::to_vector(
                          h_multiple.handles |
                          view::transform([&astor_name](auto&& h) {
                            return out_formatter(
                                h, storage::prop<astor_name.str>(h));
                          }));
                    },
                    py::return_value_policy::reference)
                .def(fmt::format("set_{}", astor_name.str.value).c_str(),
                     [&astor_name](
                         handles_wrap<std::vector<handle_t>>& h_multiple,
                         target_type val) {
                       ranges::for_each(h_multiple.handles, [&astor_name,
                                                             &val](auto&& h) {
                         in_formatter(h, storage::prop<astor_name.str>(h)) =
                             in_formatter(h, val);
                       });
                     })
                .def(fmt::format("multiple_set_{}", astor_name.str.value)
                         .c_str(),
                     [&astor_name](
                         handles_wrap<std::vector<handle_t>>& h_multiple,
                         std::vector<target_type> vals) {
                       ranges::for_each(
                           view::zip(vals, h_multiple.handles),
                           [&astor_name](auto&& pair) {
                             auto& val = std::get<0>(pair);
                             auto& h = std::get<1>(pair);
                             in_formatter(h,
                                          storage::prop<astor_name.str>(h)) =
                                 in_formatter(h, val);
                           });
                     });
          });
    }

    if constexpr (!match::without_attributes_bindings<
                      typename handle_t::type>) {
      mp::fold_left(
          pattern::attributes(typename handle_t::type{}),  // all attributes
          std::ref(pyhandle[1_c]),                         // initial state
          []<match::attribute A>(py::class_<handle_t> dc, A a) {
            using attr_value_t = std::decay_t<decltype(out_formatter(
                handle_t{},
                storage::get<A>(
                    std::declval<siconos::python::disks::idata_t&>(), 0,
                    handle_t{})))>;
            return dc
                .def(
                    fmt::format("{}", pattern::attribute_name(a)).c_str(),
                    [](handle_t& h) -> attr_value_t {
                      return out_formatter(h, storage::get<A>(h.data(), h));
                    },
                    py::return_value_policy::reference)
                .def(
                    fmt::format("{}_at_step", pattern::attribute_name(a))
                        .c_str(),
                    [](handle_t& h, indice_t step) -> attr_value_t {
                      return out_formatter(
                          h, storage::get<A>(h.data(), step, h));
                    },
                    py::return_value_policy::reference)

                .def(
                    fmt::format("set_{}", pattern::attribute_name(a)).c_str(),
                    [](handle_t& h, attr_value_t v) {
                      in_formatter(h, storage::get<A>(h.data(), h)) =
                          in_formatter(h, v);
                    })
                .def(fmt::format("set_{}_at_step", pattern::attribute_name(a))
                         .c_str(),
                     [](handle_t& h, attr_value_t v, indice_t step) {
                       in_formatter(h, storage::get<A>(h.data(), step, h)) =
                           in_formatter(h, v);
                     });
          });
    }

    if constexpr (match::batch_capable<item_t>) {
      mp::fold_left(
          pattern::attributes(typename handle_t::type{}),  // all attributes
          std::ref(pyhandle[3_c]),                         // initial state
          []<match::attribute A>(
              py::class_<handles_wrap<std::vector<handle_t>>> dc, A a) {
            using attr_value_t = std::decay_t<decltype(out_formatter(
                handle_t{},
                storage::get<A>(
                    std::declval<siconos::python::disks::idata_t&>(), 0,
                    handle_t{})))>;
            return dc
                .def("get",
                     [](handles_wrap<std::vector<handle_t>>& h_multiple) {
                       return h_multiple.handles;
                     })
                .def(
                    fmt::format("multipe_{}", pattern::attribute_name(a))
                        .c_str(),
                    [](handles_wrap<std::vector<handle_t>>& h_multiple) {
                      return ranges::to_vector(
                          h_multiple.handles | view::transform([](auto&& h) {
                            return out_formatter(
                                h, storage::get<A>(h.data(), h));
                          }));
                    },
                    py::return_value_policy::reference)
                .def(
                    fmt::format("set_{}", pattern::attribute_name(a)).c_str(),
                    [](handles_wrap<std::vector<handle_t>>& h_multiple,
                       attr_value_t val) {
                      ranges::for_each(h_multiple.handles, [&val](auto&& h) {
                        in_formatter(h, storage::get<A>(h.data(), h)) =
                            in_formatter(h, val);
                      });
                    })
                .def(
                    fmt::format("multiple_set_{}", pattern::attribute_name(a))
                        .c_str(),
                    [](handles_wrap<std::vector<handle_t>>& h_multiple,
                       std::vector<attr_value_t> vs) {
                      ranges::for_each(
                          view::zip(vs, h_multiple.handles), [](auto&& pair) {
                            auto& v = std::get<0>(pair);
                            auto& h = std::get<1>(pair);
                            in_formatter(h, storage::get<A>(h.data(), h)) =
                                in_formatter(h, v);
                          });
                    });
          });
    }

    //    using item_t = typename handle_t::type;
    mp::fold_left(storage::methods(pyhandle[2_c]),  // all methods
                  std::ref(pyhandle[1_c]),          // initial state
                  []<typename M>(py::class_<handle_t> dc, M m) {
                    return dc.def(pattern::method_name(m),
                                  pattern::method_def(m),
                                  py::return_value_policy::reference);
                  });

    auto item_name = storage::bind_name<item_t, disks_properties_t>();

    disks.def(
        fmt::format("add_{}", item_name).c_str(),
        [](siconos::python::disks::data_t& data) {
          return siconos::storage::add<item_t>(data());
        },
        py::return_value_policy::reference, R"pbdoc(
        Add
    )pbdoc");

    if constexpr (match::batch_capable<item_t>) {
      disks.def(
          fmt::format("multiple_add_{}", item_name).c_str(),
          [](siconos::python::disks::data_t& data, indice_t count) {
            auto r = siconos::storage::add<item_t>(data(), count);
            return handles_wrap<std::vector<handle_t>>{
                ranges::to<std::vector<handle_t>>(r)};
          },
          R"pbdoc(
        Add
    )pbdoc");
    }
  });

  disks.attr("__version__") = "dev";
  m.attr("__version__") = "dev";
}
