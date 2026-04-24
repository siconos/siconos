#pragma once

#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/utils/range.hpp"

namespace siconos::storage {

using namespace pattern;

/**
 * @brief Forward declaration of the handle template class
 *
 * @tparam B Base class template for the handle implementation
 * @tparam T Item type that this handle manages
 * @tparam R Index value type
 * @tparam D Data storage type
 */
template <template <typename... Args> typename B, match::item T, typename R,
          typename D>
struct handle;

/**
 * @brief Index class for uniquely identifying items in storage
 *
 * @tparam T The item type this index refers to
 * @tparam R The underlying value type for the index
 */
template <match::item T, typename R>
struct index {
  /// Tag type for index concept
  using index_t = void;
  /// Tag type for handle concept (deprecated)
  using handle_t = void;
  /// The item type this index refers to
  using type = T;
  /// The value type stored in the index
  using value_t = R;

  /// The actual index value
  value_t _value = {};

  /// @brief Get mutable reference to the index value
  /// @return Reference to the stored value
  value_t& value() noexcept { return _value; };

  /// @brief Get const reference to the index value
  /// @return Const reference to the stored value
  const value_t& value() const noexcept { return _value; };

  /// @brief Set the index value
  /// @param new_value The new value to assign
  void set_value(value_t new_value) noexcept { _value = new_value; };

  /**
   * @brief Construct index from a reference value
   * @param ref Initial index value
   */
  index(R ref) : _value{ref} {}

  /// @brief Default constructor
  index() : _value{} {};

  /**
   * @brief Constructor from a handle
   * @tparam B Base class template of the handle
   * @tparam D Data storage type
   * @param h Handle to copy index from
   */
  template <template <typename...> typename B, typename D>
  index(const handle<B, T, R, D>& h) : _value(h.index().value())
  {
  }

  /**
   * @brief Assignment operator from a handle
   * @tparam B Base class template of the handle
   * @tparam D Data storage type
   * @param h Handle to copy index from
   * @return Reference to this index
   */
  template <template <typename...> typename B, typename D>
  index<T, R>& operator=(const handle<B, T, R, D>& h)
  {
    _value = h.index().value();
    return *this;
  }

  /// @brief Get the item type associated with this index
  constexpr auto item_type() { return type{}; };

  /// @brief Three-way comparison operator
  friend auto operator<=>(const index<T, R>&, const index<T, R>&) = default;
};

/**
 * @brief Base handle implementation that stores data and index by value
 *
 * @tparam T Item type
 * @tparam R Index value type
 * @tparam D Data storage type
 */
template <typename T, typename R, typename D>
struct handle_base {
  /// Reference to the data storage
  D& _data;
  /// Index for the managed item
  index<T, R> _index;
};

/**
 * @brief Handle implementation that stores data by reference and index by
 * reference
 *
 * @tparam T Item type
 * @tparam R Index value type
 * @tparam D Data storage type
 */
template <typename T, typename R, typename D>
struct handle_ref {
  /// Reference to the data storage
  D& _data;
  /// Reference to the index
  index<T, R>& _index;
};

/**
 * @brief Main handle class for accessing items in storage
 *
 * @tparam B Base class template (either handle_base or handle_ref)
 * @tparam T Item type this handle manages
 * @tparam R Index value type
 * @tparam D Data storage type
 */
template <template <typename... Args> typename B, match::item T, typename R,
          typename D>
struct handle : B<T, R, D>, T::template interface<handle<B, T, R, D>> {
  /// Base type alias
  using base_t = B<T, R, D>;
  /// Index type for this handle
  using index_t = index<T, R>;
  /// Item type
  using type = typename index_t::type;
  /// Tag for handle concept (deprecated)
  using handle_t = void;
  /// Tag for full handle concept (deprecated)
  using full_handle_t = void;
  /// Info type extracted from data storage
  using info_t = get_info_t<D>;
  /// Data storage type
  using data_t = D;
  /// Index type from environment
  using indice = typename info_t::template env<T>::indice;
  /// Type of attached storages for this item
  using attached_storages_t = decltype(mp::filter(
      typename info_t::all_properties_t{}, mp::is_a_model<[]<typename X>() {
        return match::attached_storage<X, T>;
      }>));

  /// @brief Get the item type managed by this handle
  auto item_type() { return type{}; }

  /// @brief Get the index of this handle
  constexpr decltype(auto) index() { return base_t::_index; }

  /// @brief Get the const index of this handle
  constexpr decltype(auto) index() const { return base_t::_index; }

  /// @brief Get reference to the underlying data storage
  constexpr decltype(auto) data() { return base_t::_data; };

  /// @brief Get all attributes for this item type
  constexpr decltype(auto) attributes()
  {
    return attributes(typename index_t::type{});
  };

  /**
   * @brief Access an attached storage property by tag
   * @tparam A Tag type to identify the attached storage
   * @param step Memory step index (for time-dependent storage)
   * @return Reference to the attached storage value
   */
  template <typename A>
  constexpr decltype(auto) property(A, indice step = 0)
  {
    using item_t = T;
    // Filter properties to find attached storage matching the tag
    constexpr auto tpl = mp::filter(
        typename info_t::all_properties_t{}, mp::is_a_model<[]<typename X>() {
          return (match::attached_storage<X, item_t> && (match::tag<X, A>));
        }>);

    // Ensure we found at least one matching attached storage
    static_assert(mp::size(tpl) >= mp::size_c<1>,
                  "attached storage not found");

    using attached_storage_t = std::decay_t<decltype(tpl[0_c])>;
    return memory(step, mp::get<attached_storage_t>(
                            data().store()))[this->index().value()];
  }

  /**
   * @brief Access an attached storage property by compile-time string
   * @tparam S Compile-time string literal
   * @param step Memory step index
   * @return Reference to the attached storage value
   */
  template <string_literal S>
  constexpr decltype(auto) property(indice step = 0)
  {
    return property(symbol<S>{}, step);
  }

  /**
   * @brief Construct handle from rvalue data and rvalue reference
   * @param data Data storage
   * @param ref Index value
   */
  explicit handle(D&& data, R&& ref) : base_t{data, ref} {};

  /**
   * @brief Construct handle from lvalue data and lvalue reference
   * @param data Data storage
   * @param ref Index value
   */
  explicit handle(D& data, R& ref) : base_t{data, ref} {};

  /**
   * @brief Construct handle from lvalue data and const reference
   * @param data Data storage
   * @param ref Index value
   */
  explicit handle(D& data, const R& ref) : base_t{data, ref} {};

  /**
   * @brief Construct handle from data and existing index
   * @param data Data storage
   * @param indx Index object
   */
  explicit handle(D& data, index_t& indx) : base_t{data, indx} {};

  /**
   * @brief Construct handle from data and rvalue index
   * @param data Data storage
   * @param indx Index object (moved)
   */
  explicit handle(D& data, index_t&& indx) : base_t{data, indx} {};

  /// @brief Default constructor
  handle() : base_t{} {};

  /// @brief Copy constructor from handle reference
  handle(const handle& h) : base_t{h._data, h._index} {};

  /// @brief Move constructor
  handle(handle&&) = default;

  /**
   * @brief Assignment operator from compatible handle type
   * @tparam OtherBase Base template of the source handle
   * @param other Source handle to copy index from
   * @return Reference to this handle
   */
  template <template <typename...> typename OtherBase>
  handle& operator=(const handle<OtherBase, T, R, D>& other)
  {
    this->_index = other.index();
    return *this;
  }

  /**
   * @brief Assignment operator from same handle type
   * @param h Source handle
   * @return Reference to this handle
   */
  handle& operator=(const handle& h)
  {
    this->_index = h.index();
    return *this;
  };

  // handle operator=(handle&& h) { return handle(h); };
  /// @brief Three-way comparison operator
  friend auto operator<=>(const handle<B, T, R, D>&,
                          const handle<B, T, R, D>&) = default;
};

/**
 * @brief Create a full handle for an item
 * @tparam T Item type
 * @param data Data storage
 * @param ref Index value
 * @return Handle object
 */
template <typename T>
auto make_full_handle(auto& data, const auto& ref)
{
  using data_t = std::decay_t<decltype(data)>;
  using info_t = get_info_t<data_t>;
  using indice = typename info_t::template env<T>::indice;

  return handle<handle_base, T, indice, data_t>{data, ref};
}

/**
 * @brief Create a handle from another handle
 * @tparam T Item type
 * @param h Source handle
 * @return New handle referencing the same data
 */
template <match::item T>
auto make_handle(auto& h)
{
  return make_full_handle<T>(h.data(), h.index().value());
}

/**
 * @brief Create a handle from data and index
 * @tparam T Item type
 * @tparam R Index value type
 * @tparam D Data storage type
 * @param data Data storage
 * @param indx Index object
 * @return Handle object
 */
template <typename T, typename R, typename D>
auto make_handle(D&& data, index<T, R>&& indx)
{
  return handle<handle_base, T, R, D>{data, indx};
}

template <typename T, typename R, typename D>
auto make_handle(D& data, index<T, R>& indx)
{
  return handle<handle_base, T, R, D>{data, indx};
}



/**
 * @brief Create a handle from data and rvalue index
 * @tparam T Item type
 * @tparam R Index value type
 * @tparam D Data storage type
 * @param data Data storage
 * @param indx Index object (moved)
 * @return Handle object
 */
template <typename T, typename R, typename D>
auto make_handle(D& data, index<T, R>&& indx)
{
  return handle<handle_base, T, R, D>{data, static_cast<index<T, R>&&>(indx)};
}

/**
 * @brief Create a reference handle (stores index by reference)
 * @tparam T Item type
 * @tparam R Index value type
 * @tparam D Data storage type
 * @param data Data storage
 * @param indx Index object
 * @return Reference handle object
 */
template <typename T, typename R, typename D>
auto make_ref_handle(D& data, index<T, R>& indx)
{
  return handle<handle_ref, T, R, D>{data, indx};
}

/**
 * @brief Create a reference handle from rvalue index
 * @tparam T Item type
 * @tparam R Index value type
 * @tparam D Data storage type
 * @param data Data storage
 * @param indx Index object (moved)
 * @return Reference handle object
 */
template <typename T, typename R, typename D>
auto make_ref_handle(D& data, index<T, R>&& indx)
{
  return handle<handle_ref, T, R, D>{data, static_cast<index<T, R>&&>(indx)};
}

/**
 * @brief Generator for creating ranges of handles to all items of a given
 * type
 * @tparam I Item type
 * @param data Data storage
 * @param step Memory step for time-dependent data
 * @return Range of handles
 */
template <match::item I>
static auto handles =
    [](auto& data, std::size_t step = 0) constexpr -> decltype(auto) {
  using info_t = get_info_t<decltype(data)>;
  using env = typename info_t::template env<I>;
  using indice = typename env::indice;

  using attributes_t = std::decay_t<decltype(attributes(I{}))>;
  // Get the number of items from the first attribute's storage
  // This assumes at least one attribute exists, which may not hold for some
  // relation types
  indice num = std::size(attr_values<nth_t<0, attributes_t>>(data, step));
  return view::iota((indice)0, num) | view::transform([&data](indice i) {
           return handle<handle_base, I, indice,
                         std::decay_t<decltype(data)>>(data,
                                                       index<I, indice>(i));
         });
};

// Concept check for __init__ member function presence
static constexpr auto has_init =
    mp::is_valid([](auto&& x) -> decltype(x.__init__()) {});

// Concept check for __del__ member function presence
static constexpr auto has_del =
    mp::is_valid([](auto&& x) -> decltype(x.__del__()) {});

/**
 * @brief Concept for checking if a handle's item type derives from B
 * @tparam B Base type to check against
 */
template <typename B>
static constexpr auto handle_derive_from =
    mp::is_a_model<[]<typename T>() consteval {
      return std::derived_from<typename T::type, B>;
    }>;

/**
 * @brief Concept for checking if a handle's item type does not derive from B
 * @tparam B Base type to check against
 */
template <typename B>
static constexpr auto not_handle_derive_from =
    mp::is_a_model<[]<typename T>() consteval {
      return std::derived_from<typename T::type, B>;
    }>;

}  // namespace siconos::storage
