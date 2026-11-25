/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*! \file SpaceFilter_impl.hpp
 *  \brief implementation details for moving plans
 */

#ifndef SpaceFilter_impl_hpp
#define SpaceFilter_impl_hpp
#include <memory>

#include "SiconosSerialization.hpp"

namespace siconos::modeling {
class DynamicalSystem;
}  // namespace siconos::modeling

// Internal stuff. Not supposed to be used/viewed in users'interface.
namespace siconos::collision::native::internal {

class Hashed : public std::enable_shared_from_this<Hashed> {
 protected:
  ACCEPT_SERIALIZATION(Hashed);

  Hashed() = default;
  Hashed(const Hashed&) = delete;
  Hashed(Hashed&&) = delete;
  Hashed operator=(const Hashed&) = delete;
  Hashed operator=(Hashed&&) = delete;

 public:
  std::shared_ptr<siconos::modeling::DynamicalSystem> body{nullptr};
  int i{0};
  int j{0};
  int k{0};

  Hashed(std::shared_ptr<siconos::modeling::DynamicalSystem> body, int i, int j, int k = 0)
      : body(body), i(i), j(j), k(k) {};

  Hashed(int i, int j, int k = 0) : i(i), j(j), k(k) {};

  ~Hashed() noexcept = default;
};

bool operator==(std::shared_ptr<Hashed> const& a, std::shared_ptr<Hashed> const& b);

std::size_t hash_value(std::shared_ptr<Hashed> const& h);

}  // namespace siconos::collision::native::internal
#endif
