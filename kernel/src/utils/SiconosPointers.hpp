/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
