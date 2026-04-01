#pragma once

#include "siconos/storage/pattern/base.hpp"

namespace siconos::storage::some {

using namespace siconos::storage::pattern;

template <typename... Args>
struct attribute {
  using args = pattern::gather<Args...>;
  using attribute_t = void;
};

struct property {
  using property_t = void;
};

struct attached_storage : property {};

// not here
struct time_invariant : property {
  using time_invariant_t = void;
};

struct boolean : attribute<> {};

struct scalar : attribute<> {};

struct indice : attribute<> {};

struct integer : attribute<> {};

struct undefined_indice_parameter : attribute<> {};

template <string_literal S>
struct indice_parameter : undefined_indice_parameter {
  static constexpr auto name = S;
};

struct undefined_type_parameter : attribute<> {};

template <string_literal S>
struct type_parameter : undefined_type_parameter {
  static constexpr auto name = S;
};

template <string_literal S>
struct item : type_parameter<S> {};

struct undefined_indice_value : attribute<> {};
template <auto I>
struct indice_value : undefined_indice_value {
  static constexpr auto value = I;
};

struct undefined_vdescriptor : attribute<> {};
struct undefined_edescriptor : attribute<> {};

template <typename T>
struct vdescriptor : undefined_vdescriptor {
  using vdescriptor_t = void;
  static constexpr bool descriptor = true;
  using type = T;
};

template <typename T>
struct edescriptor : undefined_edescriptor {
  using edescriptor_t = void;
  static constexpr bool descriptor = true;
  using type = T;
};

struct unbounded_storage : attribute<> {};
struct bounded_storage : attribute<> {};
struct sparse_storage : attribute<> {};

struct undefined_unbounded_collection : unbounded_storage {};
struct undefined_bounded_collection : bounded_storage {};

struct undefined_matrix : bounded_storage {};
struct undefined_diagonal_matrix : undefined_matrix {};
struct undefined_vector : bounded_storage {};
struct undefined_array : bounded_storage {};

struct undefined_unbounded_matrix : unbounded_storage {};
struct undefined_unbounded_vector : unbounded_storage {};
struct undefined_unbounded_array : unbounded_storage {};

struct undefined_sparse_matrix : sparse_storage {};

template <typename... Sizes>
struct with_sizes {
  using sizes = pattern::gather<Sizes...>;
};

template <typename Type>
struct with_type {
  using type = Type;
};

template <typename... Types>
struct with_types {
  using types = pattern::gather<Types...>;
};

template <typename Type, typename N, typename M>
struct matrix : undefined_matrix, with_sizes<N, M>, with_type<Type> {};

template <typename Mat>
struct transposed_matrix
    : undefined_matrix,
      with_sizes<nth_t<1, typename Mat::sizes>, nth_t<0, typename Mat::sizes>>,
      with_type<typename Mat::type> {};

template <typename Type = some::scalar>
struct unbounded_matrix : undefined_unbounded_matrix, with_type<Type> {};

template <typename Type, typename N>
struct unbounded_col_matrix : undefined_unbounded_matrix, with_type<Type>, with_sizes<N> {};

template <typename Type = some::scalar>
struct unbounded_vector : undefined_unbounded_vector, with_type<Type> {};

template <typename Type, typename M>
struct diagonal_matrix : undefined_diagonal_matrix,
                         with_sizes<nth_t<0, typename M::sizes>>,
                         with_type<Type> {};

template <typename Type>
struct unbounded_diagonal_matrix : unbounded_storage,
                                   undefined_diagonal_matrix,
                                   with_type<Type> {};

template <typename Type>
struct sparse_matrix : undefined_sparse_matrix, with_type<Type> {};

struct structure_for_assembled_data {};

template <typename Type = some::scalar>
struct assembled_matrix : unbounded_matrix<Type>, structure_for_assembled_data {};

template <typename Type = some::scalar>
struct assembled_diagonal_matrix : unbounded_diagonal_matrix<Type>,
                                   structure_for_assembled_data {};

template <typename Type = some::scalar>
struct assembled_vector : unbounded_vector<Type>, structure_for_assembled_data {};

template <typename Type, typename N>
struct vector : undefined_vector, with_sizes<N>, with_type<Type> {};

template <typename Type, typename N>
struct array : undefined_array, with_sizes<N>, with_type<Type> {};

template <typename Type>
struct unbounded_array : undefined_unbounded_array, with_type<Type> {};

struct undefined_map : attribute<> {};
template <typename Key, typename Value>
struct map : undefined_map, with_types<Key, Value> {};

struct undefined_graph : attribute<> {};

template <typename Edge, typename Vertice>
struct graph : undefined_graph, with_types<Edge, Vertice> {};

template <typename Type>
struct unbounded_collection : undefined_unbounded_collection, with_type<Type> {};

template <typename Type, typename N>
struct bounded_collection : undefined_bounded_collection, with_sizes<N>, with_type<Type> {};

// xx should be elsewhere
template <match::item T>
struct item_ref : attribute<>, with_type<T> {};

struct undefined_polymorphic_type {};
template <typename... Ts>
struct polymorph : undefined_polymorphic_type, with_types<Ts...> {
  using polymorphic = void;
};

template <typename... Ts>
struct polymorphic_attribute : attribute<>, with_type<gather<Ts...>>, polymorph<Ts...> {};

struct given_type {};

template <typename T>
struct specific : given_type, attribute<> {
  using xtype = T;
  using type_t = void;
};

struct given_definition : attribute<> {};

template <string_literal S>
struct definition : given_definition, symbol<S> {};

struct undefined_tuple {};

template <typename... Ts>
struct tuple : attribute<>, undefined_tuple, with_types<Ts...> {};
}  // namespace siconos::storage::some
