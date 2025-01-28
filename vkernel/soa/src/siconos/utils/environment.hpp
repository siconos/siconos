#pragma once

#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>
#include <tuple>
#include <variant>
#include <vector>
#include <map>

#include "siconos/algebra/eigen.hpp"
#include "siconos/algebra/linear_algebra.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos {
template <typename Params>
struct standard_environment {
  using params = Params;

  using boolean = uint8_t; // not bool => cf CompactNSearch sort.
  using scalar = double;
  using indice = std::size_t;
  using integer = std::int64_t;

  template <typename K, typename V>
  using map = std::map<K, V>;

  template <typename... Ts>
  using tuple = std::tuple<Ts...>;

  template <typename S>
  using unbounded_collection = std::vector<S>;

  template <typename S, std::size_t N>
  using bounded_collection = std::array<S, N>;

  template <typename T, indice M>
  using array = std::array<T, M>;

  template <typename T, indice M>
  using vector = algebra::vector<T, M>;

  template <typename T, indice M, indice N>
  using matrix = algebra::matrix<T, M, N>;

  template <typename T>
  using unbounded_vector = algebra::unbounded_vector<T>;

  template <typename T>
  using unbounded_matrix = algebra::unbounded_matrix<T>;

  template <typename T, indice M>
  using unbounded_col_matrix = algebra::unbounded_col_matrix<T, M>;

  template <typename T, indice N>
  using diagonal_matrix = algebra::diagonal_matrix<T, N>;

  template <typename T>
  using assembled_matrix = algebra::mat<T>;

  template <typename T>
  using assembled_vector = algebra::vec<T>;

  template <typename T>
  using assembled_diagonal_matrix = algebra::diag_mat<T>;

  template <typename T>
  using item_ref = storage::index<T, indice>;

  template <typename... Ts>
  using variant = std::variant<Ts...>;

  template <typename T>
  using default_storage = std::array<T, 1>;

};
}  // namespace siconos
