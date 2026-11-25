// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef SiconosPointers_hpp
#define SiconosPointers_hpp

/*! \file SiconosPointers.hpp
 */

#include <memory>  // Def of std::shared_ptr

namespace siconos::pointers {

namespace internal {

/** Using a shared_ptr to hold a pointer to a statically allocated
 object
 use create<type>SPtr(<type> &x)
 cf http://www.boost.org/doc/
*/
struct nullDeleter {
  void operator()(void const*) const {}
};

}  // namespace internal

/** create a shared_ptr from a statically allocated object, without copy.
    This should be avoided as much as possible ...
 */
template <typename T>
std::shared_ptr<T> createSPtr(T& input) {
  std::shared_ptr<T> px(&input, siconos::pointers::internal::nullDeleter());
  return px;
}
}  // namespace siconos::pointers
#endif /* SiconosPointers_hpp */
