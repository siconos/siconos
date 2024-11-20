#pragma once

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/utils/range.hpp"

namespace siconos::storage {

using namespace pattern;

template <match::item T, typename R>
struct index {
  using index_t = void;
  using handle_t = void;  // remove this
  using type = T;
  using value_t = R;

  value_t _value = {};

  value_t get() const noexcept { return _value; };

  explicit index(R ref) : _value{ref} {}

  index() : _value{} {};

  constexpr auto item_type() { return type{}; };

  friend auto operator<=>(const index<T, R>&, const index<T, R>&) = default;
};

template <match::item T, typename R, typename D>
struct handle : index<T, R>, T::template interface<handle<T, R, D>> {
  using base_index_t = index<T, R>;
  using full_handle_t = void;
  using info_t = get_info_t<D>;
  using data_t = D;
  using indice = typename info_t::env::indice;
  using attached_storages_t =
      decltype(ground::filter(typename info_t::all_properties_t{},
                              ground::is_a_model<[]<typename X>() {
                                return match::attached_storage<X, T>;
                              }>));

  data_t& _data;

  constexpr decltype(auto) index_cast()
  {
    return static_cast<index<T, R>>(*this);
  }

  constexpr decltype(auto) attributes()
  {
    return attributes(typename index<T, R>::type{});
  };

  template <typename A>
  constexpr decltype(auto) property(A, indice step = 0)
  {
    using item_t = T;
    constexpr auto tpl = ground::filter(
        typename info_t::all_properties_t{},
        ground::is_a_model<[]<typename X>() {
          return (match::attached_storage<X, item_t> && (match::tag<X, A>));
        }>);

    static_assert(ground::size(tpl) >= ground::size_c<1>,
                  "attached storage not found");

    using attached_storage_t = std::decay_t<decltype(tpl[0_c])>;
    return memory(
        step, ground::get<attached_storage_t>(data().store()))[this->get()];
  }

  // not convenient, it needs to specify template keyword:
  // some_handle.template property<S>() = ...
  template <string_literal S>
  constexpr decltype(auto) property(indice step = 0)
  {
    return property(symbol<S>{}, step);
  }

  decltype(auto) data() { return _data; };

  explicit handle(D& data, R& ref) : index<T, R>{ref}, _data{data} {};

  explicit handle(D& data, index<T, R>& ha) : index<T, R>{ha}, _data{data} {};

  explicit handle(D& data, index<T, R>&& ha)
      : index<T, R>{ha}, _data{data} {};

  handle() : index<T, R>{}, _data{} {};

  handle(const handle& h) : index<T, R>(h.get()), _data(h._data){};
  handle(handle&&) = default;
  handle operator=(const handle& h) { return handle(h); };
  handle operator=(handle&& h) { return handle(h); };
  ;

  friend auto operator<=>(const handle<T, R, D>&,
                          const handle<T, R, D>&) = default;
};

template <match::item T>
static decltype(auto) make_half_handle(auto indx)
{
  using indice = std::decay_t<decltype(indx)>;
  return index<T, indice>{indx};
}

template <typename T>
auto make_full_handle(auto& data, const auto& indx)
{
  using data_t = std::decay_t<decltype(data)>;
  using info_t = get_info_t<data_t>;
  using indice = typename info_t::env::indice;

  indice index = indx;
  return handle<T, indice, data_t>{data, index};
}

template <match::item T>
decltype(auto) make_handle(auto& h)
{
  return make_full_handle<T>(h.data(), h.get());
}

template <match::item I>
static auto handles =
    [](auto& data, std::size_t step = 0) constexpr -> decltype(auto) {
  using info_t = get_info_t<decltype(data)>;
  using env = typename info_t::env;
  using indice = typename env::indice;

  using attributes_t = std::decay_t<decltype(attributes(I{}))>;
  // need at least one attributes
  // empty items are not supposed to exist
  indice num = std::size(attr_values<nth_t<0, attributes_t>>(data, step));
  return view::iota((indice)0, num) | view::transform([&data](indice i) {
           return handle<I, indice, std::decay_t<decltype(data)>>(
               data, index<I, indice>(i));
         });
};
}  // namespace siconos::storage
