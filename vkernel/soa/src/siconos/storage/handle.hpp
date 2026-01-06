#pragma once

#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/utils/range.hpp"

namespace siconos::storage {

using namespace pattern;

template <template <typename... Args> typename B, match::item T, typename R,
          typename D>
struct handle;

template <match::item T, typename R>
struct index {
  using index_t = void;
  using handle_t = void;  // remove this
  using type = T;
  using value_t = R;

  value_t _value = {};
  value_t& value() noexcept { return _value; };
  const value_t& value() const noexcept { return _value; };
  void set_value(value_t new_value) noexcept { _value = new_value; };

  index(R ref) : _value{ref} {}
  index() : _value{} {};
  // Constructor from handle
  template <template <typename...> typename B, typename D>
  index(const handle<B, T, R, D>& h) : _value(h.index().value())
  {
  }

  // Assignment from handle
  template <template <typename...> typename B, typename D>
  index<T, R>& operator=(const handle<B, T, R, D>& h)
  {
    _value = h.index().value();
    return *this;
  }

  constexpr auto item_type() { return type{}; };

  friend auto operator<=>(const index<T, R>&, const index<T, R>&) = default;
};

template <typename T, typename R, typename D>
struct handle_base {
  D& _data;
  index<T, R> _index;
};

template <typename T, typename R, typename D>
struct handle_ref {
  D& _data;
  index<T, R>& _index;
};

template <template <typename... Args> typename B, match::item T, typename R,
          typename D>
struct handle : B<T, R, D>, T::template interface<handle<B, T, R, D>> {
  using base_t = B<T, R, D>;
  using index_t = index<T, R>;
  using type = typename index_t::type;
  using handle_t = void;
  using full_handle_t = void;
  using info_t = get_info_t<D>;
  using data_t = D;
  using indice = typename info_t::env::indice;
  using attached_storages_t = decltype(mp::filter(
      typename info_t::all_properties_t{}, mp::is_a_model<[]<typename X>() {
        return match::attached_storage<X, T>;
      }>));

  auto item_type() { return type{}; }

  constexpr decltype(auto) index() { return base_t::_index; }

  constexpr decltype(auto) index() const { return base_t::_index; }

  decltype(auto) data() { return base_t::_data; };

  constexpr decltype(auto) attributes()
  {
    return attributes(typename index_t::type{});
  };

  template <typename A>
  constexpr decltype(auto) property(A, indice step = 0)
  {
    using item_t = T;
    constexpr auto tpl = mp::filter(
        typename info_t::all_properties_t{}, mp::is_a_model<[]<typename X>() {
          return (match::attached_storage<X, item_t> && (match::tag<X, A>));
        }>);

    static_assert(mp::size(tpl) >= mp::size_c<1>,
                  "attached storage not found");

    using attached_storage_t = std::decay_t<decltype(tpl[0_c])>;
    return memory(step, mp::get<attached_storage_t>(
                            data().store()))[this->index().value()];
  }

  // not convenient, it needs to specify template keyword:
  // some_handle.template property<S>() = ...
  template <string_literal S>
  constexpr decltype(auto) property(indice step = 0)
  {
    return property(symbol<S>{}, step);
  }

  explicit handle(D& data, R& ref) : base_t{data, ref} {};

  explicit handle(D& data, const R& ref) : base_t{data, ref} {};

  explicit handle(D& data, index_t& indx) : base_t{data, indx} {};

  explicit handle(D& data, index_t&& indx) : base_t{data, indx} {};

  handle() : base_t{} {};

  handle(handle& h) : base_t{h._data, h._index} {};
  handle(handle&&) = default;

  template <template <typename...> typename OtherBase>
  handle& operator=(const handle<OtherBase, T, R, D>& other)
  {
    this->_index = other.index();
    return *this;
  }

  // Fixed existing assignment operator
  handle& operator=(const handle& h)
  {
    this->_index = h.index();
    return *this;
  };

  // handle operator=(handle&& h) { return handle(h); };
  ;
  friend auto operator<=>(const handle<B, T, R, D>&,
                          const handle<B, T, R, D>&) = default;
};

template <typename T>
auto make_full_handle(auto& data, const auto& ref)
{
  using data_t = std::decay_t<decltype(data)>;
  using info_t = get_info_t<data_t>;
  using indice = typename info_t::env::indice;

  return handle<handle_base, T, indice, data_t>{data, ref};
}

template <match::item T>
auto make_handle(auto& h)
{
  return make_full_handle<T>(h.data(), h.index().value());
}

template <typename T, typename R, typename D>
auto make_handle(D& data, index<T, R>& indx)
{
  return handle<handle_base, T, R, D>{data, indx};
}

template <typename T, typename R, typename D>
auto make_handle(D& data, index<T, R>&& indx)
{
  return handle<handle_base, T, R, D>{data, static_cast<index<T, R>&&>(indx)};
}

template <typename T, typename R, typename D>
auto make_ref_handle(D& data, index<T, R>& indx)
{
  return handle<handle_ref, T, R, D>{data, indx};
}

template <typename T, typename R, typename D>
auto make_ref_handle(D& data, index<T, R>&& indx)
{
  return handle<handle_ref, T, R, D>{data, static_cast<index<T, R>&&>(indx)};
}

template <match::item I>
static auto handles =
    [](auto& data, std::size_t step = 0) constexpr -> decltype(auto) {
  using info_t = get_info_t<decltype(data)>;
  using env = typename info_t::env;
  using indice = typename env::indice;

  using attributes_t = std::decay_t<decltype(attributes(I{}))>;
  // Warning! need at least one attributes!
  // This won't work for items without attributes such as some relations.
  indice num = std::size(attr_values<nth_t<0, attributes_t>>(data, step));
  return view::iota((indice)0, num) | view::transform([&data](indice i) {
           return handle<handle_base, I, indice,
                         std::decay_t<decltype(data)>>(data,
                                                       index<I, indice>(i));
         });
};

// Check if a  __init__ member function is present.
static constexpr auto has_init =
    mp::is_valid([](auto&& x) -> decltype(x.__init__()) {});

// Check if a  __del__ member function is present.
static constexpr auto has_del =
    mp::is_valid([](auto&& x) -> decltype(x.__del__()) {});

template <typename B>
static constexpr auto handle_derive_from =
    mp::is_a_model<[]<typename T>() consteval {
      return std::derived_from<typename T::type, B>;
    }>;

template <typename B>
static constexpr auto not_handle_derive_from =
    mp::is_a_model<[]<typename T>() consteval {
      return std::derived_from<typename T::type, B>;
    }>;

}  // namespace siconos::storage
