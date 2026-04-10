#pragma once

#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"

namespace siconos::storage {
using namespace pattern;

template <match::item I, typename Tag, match::attribute DataSpec>
struct attached : DataSpec, some::attached_storage {
  using item = I;
  using tag = Tag;
  using data_spec = DataSpec;
};

namespace property {
struct keep : some::property {};

struct wrapped : some::property {};

struct time_invariant : some::property {};

struct refine : some::property {};

struct bind : some::property {};

struct diagonal : refine {
  template <match::attribute A>
  using refine = some::diagonal_matrix<typename A::type, A>;
};
struct assembled_diagonal : refine {
  template <match::attribute A>
  using refine = some::assembled_diagonal_matrix<typename A::type>;
};

struct unbounded_vector : refine {
  template <match::attribute A>
  using refine = some::unbounded_vector<typename A::type>;
};

struct unbounded_matrix : refine {
  template <match::attribute A>
  using refine = some::unbounded_matrix<typename A::type>;
};

struct assembled_vector : refine {
  template <match::attribute A>
  using refine = some::assembled_vector<typename A::type>;
};

struct assembled_matrix : refine {
  template <match::attribute A>
  using refine = some::assembled_matrix<typename A::type>;
};

struct sparse_matrix : refine {
  template <match::attribute A>
  using refine = some::sparse_matrix<typename A::type>;
};

}  // namespace property

template <match::attribute A, typename T>
struct refine_with_type : A {
  using type = T;
};

template <match::property... Parts>
struct with_properties : item {
  using properties = gather<Parts...>;
};

template <match::attribute Attr, std::size_t N>
struct keep : property::keep {
  using type = Attr;
  using keep_t = void;
  //    using attribute = Attr;
  static constexpr std::size_t size = N;
};

template <match::item Item, template <typename... Ts> typename Wrapper>
struct wrapped : property::wrapped {
  using wrapped_t = void;
  template <typename... Ts>
  using wrapper = Wrapper<Ts...>;
  using type = Item;
};

template <match::attribute Attr>
struct time_invariant : property::time_invariant {
  using type = Attr;
  using time_invariant_t = void;
};

template <match::item Item, string_literal S>
struct bind : property::bind, symbol<S> {
  using item = Item;
  using bind_t = void;
};

template <match::abstract_matrix M>
struct diagonal : property::diagonal {
  using type = M;
  using diagonal_t = void;
};

template <match::abstract_matrix M>
struct assembled_diagonal : property::assembled_diagonal {
  using type = M;
  using assembled_diagonal_t = void;
};

template <typename T>
struct unbounded;

template <match::abstract_matrix M>
struct unbounded<M> : property::unbounded_matrix {
  using type = M;
  using unbounded_matrix_t = void;
};

template <match::abstract_vector V>
struct unbounded<V> : property::unbounded_vector {
  using type = V;
  using unbounded_vector_t = void;
};

template <typename T>
struct assembled;

template <match::abstract_matrix M>
struct assembled<M> : property::assembled_matrix {
  using type = M;
  using assembled_matrix_t = void;
};

template <match::abstract_vector V>
struct assembled<V> : property::assembled_vector {
  using type = V;
  using assembled_vector_t = void;
};

template <match::abstract_matrix M>
struct sparse : property::sparse_matrix {
  using type = M;
  using sparse_matrix_t = void;
};

template <match::property K>
static auto pre_map_all_properties_as =
    []<typename D>(D& data) constexpr -> auto {
  using info_t = get_info_t<D>;
  using all_properties_t = typename info_t::all_properties_t;

  return mp::filter(all_properties_t{}, mp::derive_from<K>);
};

template <match::property K>
static auto all_properties_as = []<typename D>(D& data) constexpr -> auto {
  using info_t = get_info_t<D>;
  using all_properties_t = typename info_t::all_properties_t;

  return mp::filter(all_properties_t{}, mp::derive_from<K>);
};

template <match::attribute Attr>
static auto attribute_properties = [](auto& data) constexpr -> auto {
  using info_t = get_info_t<decltype(data)>;
  using all_properties_t = typename info_t::all_properties_t;

  return filter<hold<decltype([]<typename T>(T) {
    return std::derived_from<typename T::type, Attr>;
  })>>(all_properties_t{});
};

template <match::item Item, typename properties>
static constexpr auto item_properties_from()
{
  return mp::filter(properties{}, mp::is_a_model<[]<typename T>() consteval {
                      if constexpr (match::item_property<T>) {
                        return std::derived_from<Item, typename T::item>;
                      }
                      else {
                        return false;
                      };
                    }>);
};

template <match::item Item, typename Data>
static constexpr auto item_properties(Data&& data)
{
  using info_t = get_info_t<Data>;
  using all_properties_t = typename info_t::all_properties_t;

  return item_properties_from<Item, all_properties_t>();
};

template <match::attribute Attr, match::property K>
static constexpr bool has_property(auto& data)
{
  return mp::any_of(all_properties_as<K>(data), []<match::property P>(P) {
    return std::derived_from<Attr, typename P::type>;
  });
}

template <match::item Item, match::property K, typename properties>
static constexpr bool has_property_from()
{
  return mp::any_of(
      item_properties_from<Item, properties>(),
      []<match::property P>(P) { return std::derived_from<P, K>; });
};

template <match::item Item, match::property K>
static constexpr bool has_property(auto&& data)
{
  return mp::any_of(item_properties<Item>(data), []<match::property P>(P) {
    return std::derived_from<P, K>;
  });
};

template <match::item Item, typename Properties>
static constexpr auto bind_name()
{
  return mp::find_if(item_properties_from<Item, Properties>(),
                     mp::is_a_model<[]<match::property P>() {
                       return std::derived_from<P, property::bind>;
                     }>)
      .value_or([]<bool flag = false>() {
        static_assert(flag, "no binding found!");
      })
      .str.value;
};

template <typename A, typename K, typename D>
using has_property_t = std::decay_t<decltype(has_property<A, K>(D{}))>;

static auto refine_attribute = []<match::attribute Attr, typename D>(
                                   const D& data,
                                   Attr) constexpr -> decltype(auto) {
  using refines =
      decltype(mp::filter(pre_map_all_properties_as<property::refine>(data),
                          mp::is_inside_type_parent<Attr>));

  if constexpr (mp::size(refines{}) > mp::size_c<0_c>) {
    return typename nth_t<0, refines>::template refine<Attr>{};
  }
  else {
    // return attribute as it is
    return Attr{};
  }
};

// recursive type refinement
template <typename Attr>
struct recursive_rebuild {
  using type = Attr;
};

// specialization for internal attribute types
template <typename Attr>
  requires(match::attribute_with_internal_type<Attr> &&
           match::attribute<typename Attr::type>)
struct recursive_rebuild<Attr> {
  using inner = typename recursive_rebuild<typename Attr::type>::type;
  using type = refine_with_type<Attr, inner>;
};

static constexpr auto refine_recursively_attribute =
    []<match::attribute Attr>(auto& data, Attr) {
      using data_t = std::decay_t<decltype(data)>;

      // rebuild
      using recursed = typename recursive_rebuild<Attr>::type;

      // user-defined refinements at the top level
      return refine_attribute(data_t{}, recursed{});
    };

template <typename Handle, typename Data>
using attached_storages_t =
    std::decay_t<decltype(attached_storages(Handle{}.item_type(), Data{}))>;

template <typename Item>
constexpr decltype(auto) attached_storages(Item, auto& data)
{
  using info_t = get_info_t<decltype(data)>;
  using item_t = Item;

  return mp::filter(typename info_t::all_properties_t{},
                    mp::is_a_model<[]<typename T>() {
                      return match::attached_storage<T, item_t>;
                    }>);
};

template <typename Item>
constexpr decltype(auto) all_storages(Item, auto& data)
{
  using item_t = Item;
  return mp::tuple_unique(
      concat(attributes(item_t{}), attached_storages(item_t{}, data)));
}

// use storage::attached_storages(...) instead
template <typename Item>
static constexpr auto is_attached_storage =
    mp::is_a_model<[]<typename T>() constexpr {
      return match::attached_storage<T, Item>;
    }>;

}  // namespace siconos::storage
