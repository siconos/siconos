/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2025 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*! \file StorageTools.hpp
  Tools to handle multiple types of matrix storage (e.g storage for Dense, view ...)
*/

#ifndef STORAGE_TOOLS_H
#define STORAGE_TOOLS_H

#include <variant>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::algebra {

using OwnedDenseVector = std::unique_ptr<std::vector<double>>;
using MapDenseVector = std::shared_ptr<siconos::algebra::MapVectorType>;
using DenseVectorStorage = std::variant<std::monostate, OwnedDenseVector, MapDenseVector>;

using OwnedSparse = std::unique_ptr<siconos::algebra::SiconosSparseMatrix>;
using MapSparse = std::shared_ptr<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>;
using SparseStorage = std::variant<std::monostate, OwnedSparse, MapSparse>;

enum class AccessMode { OwnedOnly, Any };

/**
 * @brief Visit and operate on a sparse matrix storage variant
 *
 * This helper function dispatches a callable `f` to the concrete type
 * currently held by a sparse matrix storage variant. The storage variant
 * can contain either:
 *   - a `std::monostate` (uninitialized),
 *   - an `OwnedSparse` (unique_ptr to a sparse matrix),
 *   - or a `MapSparse` (shared_ptr to an Eigen::Map of a sparse matrix).
 *
 * The function enforces access rules defined by the template parameter
 * `mode`, which can restrict access to owned matrices only.
 *
 * @tparam mode
 *   Access control policy:
 *   - `AccessMode::Any`: allows visiting both owned and mapped storages.
 *   - `AccessMode::OwnedOnly`: restricts access to `OwnedSparse` only.
 *
 * @tparam Variant
 *   A `std::variant` type compatible with `SparseStorage`.
 *
 * @tparam F
 *   A callable type taking a reference to a `SiconosSparseMatrix`
 *   (either mutable or const, depending on the overload).
 *
 * @param[in,out] v
 *   Variant holding the sparse matrix storage.
 *
 * @param[in] f
 *   Function or lambda to invoke with the dereferenced matrix.
 *
 * @param[in] name
 *   Optional identifier used in error messages (default: `"storage"`).
 *
 * @return
 *   The return value of `f`, if it returns something; otherwise `void`.
 *
 * @throws std::runtime_error
 *   If the storage is unset (`std::monostate`) or if the access mode
 *   forbids the current storage type.
 *
 */
template <AccessMode mode = AccessMode::Any, typename Variant, typename F>
inline decltype(auto) visitSparseStorage(Variant& v, F&& f, const char* name = "storage") {
  using ReturnType =
      decltype(f(*std::declval<typename std::add_lvalue_reference_t<
                     typename std::variant_alternative_t<1, std::decay_t<Variant>>>>()));

  if constexpr (std::is_void_v<ReturnType>) {
    std::visit(
        [&](auto& s) {
          using T = std::decay_t<decltype(s)>;

          if constexpr (std::is_same_v<T, std::monostate>) {
            throw std::runtime_error(std::string(name) + " is not set");
          } else if constexpr (mode == AccessMode::OwnedOnly) {
            if constexpr (std::is_same_v<T, OwnedSparse>)
              f(*s);
            else
              throw std::runtime_error(std::string(name) + " must contain an OwnedSparse");
          } else {  // AccessMode::Any
            f(*s);
          }
        },
        v);
  } else {
    return std::visit(
        [&](auto& s) -> ReturnType {
          using T = std::decay_t<decltype(s)>;

          if constexpr (std::is_same_v<T, std::monostate>) {
            throw std::runtime_error(std::string(name) + " is not set");
          } else if constexpr (mode == AccessMode::OwnedOnly) {
            if constexpr (std::is_same_v<T, OwnedSparse>)
              return f(*s);
            else
              throw std::runtime_error(std::string(name) + " must contain an OwnedSparse");
          } else {  // AccessMode::Any
            return f(*s);
          }
        },
        v);
  }
}

template <AccessMode mode = AccessMode::Any, typename Variant, typename F>
inline decltype(auto) visitSparseStorage(const Variant& v, F&& f,
                                         const char* name = "storage") {
  using ReturnType = decltype(f(
      *std::declval<const typename std::variant_alternative_t<1, std::decay_t<Variant>>&>()));

  if constexpr (std::is_void_v<ReturnType>) {
    std::visit(
        [&](const auto& s) {
          using T = std::decay_t<decltype(s)>;

          if constexpr (std::is_same_v<T, std::monostate>) {
            throw std::runtime_error(std::string(name) + " is not set");
          } else if constexpr (mode == AccessMode::OwnedOnly) {
            if constexpr (std::is_same_v<T, OwnedSparse>)
              f(*s);
            else
              throw std::runtime_error(std::string(name) + " must contain an OwnedSparse");
          } else {
            f(*s);
          }
        },
        v);
  } else {
    return std::visit(
        [&](const auto& s) -> ReturnType {
          using T = std::decay_t<decltype(s)>;

          if constexpr (std::is_same_v<T, std::monostate>) {
            throw std::runtime_error(std::string(name) + " is not set");
          } else if constexpr (mode == AccessMode::OwnedOnly) {
            if constexpr (std::is_same_v<T, OwnedSparse>)
              return f(*s);
            else
              throw std::runtime_error(std::string(name) + " must contain an OwnedSparse");
          } else {
            return f(*s);
          }
        },
        v);
  }
}

/**
 * @brief Unified helper to compute or update a sparse matrix via C++ or Python callbacks.
 *
 * This utility function factorizes the logic used in multiple `computeXXX`
 * methods of `LagrangianSparseDS` and similar classes.
 * It invokes either a native C++ function (`func_cpp`) or a Python callback
 * (`func_py`) depending on which is available.
 *
 * The function automatically unwraps the sparse matrix from its variant
 * storage (via `visitSparseStorage`) and applies the provided computation.
 *
 * As an example, check computeMass in LagrangianSparseDS class
 *
 * @tparam Storage
 *   Type of the sparse matrix storage (typically `SparseStorage`).
 *
 * @tparam FuncCpp
 *   Type of the C++ callback, e.g. `std::function<void(..., SiconosSparseMatrix&)>`.
 *
 * @tparam FuncPy
 *   Type of the Python callback, e.g. `std::function<SiconosSparseMatrix(...,
 * SiconosSparseMatrix&)>`.
 *
 * @tparam Args
 *   Variadic argument pack forwarded to the callbacks.
 *
 * @param[in,out] storage
 *   Variant holding the sparse matrix to be updated.
 *
 * @param[in] func_cpp
 *   Native C++ function to compute or update the sparse matrix in place.
 *
 * @param[in] func_py
 *   Python function that returns a new sparse matrix.
 *
 * @param[in] storage_name
 *   Optional name used in exception messages for clarity.
 *
 * @param[in] args
 *   Additional arguments forwarded to the computation functions
 *   (e.g., position, velocity, time).
 *
 * @return
 *   `true` if a computation was performed (either via C++ or Python),
 *   `false` if neither callback is set.
 *
 * @throws std::runtime_error
 *   If the storage is unset or incompatible with `AccessMode::OwnedOnly`.
 *
 * @details
 *   - The function ensures type safety by enforcing `AccessMode::OwnedOnly`
 *     during storage visitation.
 *   - When `func_cpp` is used, the update occurs in place.
 *   - When `func_py` is used, the storage is replaced with a new
 *     `OwnedSparse` created from the returned matrix.
 *
 */
template <typename Storage, typename FuncCpp, typename FuncPy, typename... Args>
bool computeSparse(Storage& storage, FuncCpp&& func_cpp, FuncPy&& func_py,
                   const char* storage_name, Args&&... args) {
  using siconos::algebra::AccessMode;
  using siconos::algebra::visitSparseStorage;
  if (func_cpp) {
    visitSparseStorage<AccessMode::OwnedOnly>(
        storage, [&](auto& m) { func_cpp(std::forward<Args>(args)..., m); }, storage_name);
    return true;
  } else if (func_py) {
    visitSparseStorage<AccessMode::OwnedOnly>(
        storage,
        [&](auto& m) {
          storage = std::make_unique<siconos::algebra::SiconosSparseMatrix>(
              func_py(std::forward<Args>(args)..., m));
        },
        storage_name);
    return true;
  }
  return false;
}

struct SparseSnapshot {
  const double* values_ptr;
  const int* inner_ptr;
  const int* outer_ptr;
  Eigen::Index nonzeros;
  Eigen::Index rows, cols;

  SparseSnapshot(const Eigen::SparseMatrix<double, Eigen::ColMajor>& M)
      : values_ptr(M.valuePtr()),
        inner_ptr(M.innerIndexPtr()),
        outer_ptr(M.outerIndexPtr()),
        nonzeros(M.nonZeros()),
        rows(M.rows()),
        cols(M.cols()) {}
};

inline bool isStructureStable(const SparseSnapshot& before,
                              const Eigen::SparseMatrix<double, Eigen::ColMajor>& after) {
  return before.values_ptr == after.valuePtr() && before.inner_ptr == after.innerIndexPtr() &&
         before.outer_ptr == after.outerIndexPtr() && before.nonzeros == after.nonZeros() &&
         before.rows == after.rows() && before.cols == after.cols();
}

// template <typename Variant>
// const siconos::algebra::SiconosSparseMatrix& getOwnedSparse(const Variant& storage,
//                                                             const char* name = "storage") {
//   if (!std::holds_alternative<OwnedSparse>(storage)) {
//     throw std::runtime_error(std::string(name) +
//                              " can only be accessed when it contains an OwnedSparse");
//   }

//   return std::visit(
//       [](const auto& s) -> const SiconosSparseMatrix& {
//         using T = std::decay_t<decltype(s)>;
//         if constexpr (std::is_same_v<T, OwnedSparse>) {
//           return *s;
//         } else {
//           throw std::runtime_error("Unexpected alternative in variant");
//         }
//       },
//       storage);
// }

// template <typename Variant>
// const SiconosSparseMatrix& getOwnedSparseIf(const Variant& storage,
//                                             const char* name = "storage") {
//   auto p = std::get_if<OwnedSparse>(&storage);
//   if (!p || !*p) {
//     throw std::runtime_error(std::string(name) + " must contain an OwnedSparse");
//   }
//   return **p;
// }

// template <AccessMode mode = AccessMode::OwnedOnly, typename Variant, typename F>
// decltype(auto) accessSparseStorage(const Variant& storage, F&& f,
//                                    const char* name = "storage") {
//   if constexpr (mode == AccessMode::OwnedOnly) {
//     // Access only to OwnedSparse
//     auto p = std::get_if<OwnedSparse>(&storage);
//     if (!p || !*p) {
//       throw std::runtime_error(std::string(name) + " must contain an OwnedSparse");
//     }
//     return f(**p);
//   } else {
//     // Access to both types (Owned and Map)
//     if (std::holds_alternative<std::monostate>(storage)) {
//       throw std::runtime_error(std::string(name) + " is not set");
//     }

//     if (auto p = std::get_if<OwnedSparse>(&storage); p && *p) {
//       return f(**p);
//     }
//     if (auto p = std::get_if<MapSparse>(&storage); p && *p) {
//       return f(**p);
//     }
//     throw std::runtime_error(std::string(name) + " contains unexpected alternative");
//   }
// }

}  // namespace siconos::algebra

#endif