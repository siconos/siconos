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
/*! \file RelationCompatibility.hpp
 *  \brief Helpers to keep the deprecated names of the relation API working.
 *
 *  When a virtual function of a relation is renamed (e.g. LagrangianScleronomousR::
 *  computeJacobianhOver_q -> computeH), a user class that still overrides the old name must
 *  keep being called. The pattern is:
 *
 *  \code
 *  void Base::computeH(args) {          // new name, the one called by the kernel
 *    if (!computeHLegacy_) {            // first entry: give a chance to an old override
 *      LegacyCallGuard guard(computeHLegacy_);
 *      SICONOS_IGNORE_DEPRECATED_BEGIN
 *      computeJacobianhOver_q(args);    // virtual: user override, or Base version below
 *      SICONOS_IGNORE_DEPRECATED_END
 *      return;
 *    }
 *    ...default implementation...
 *  }
 *  void Base::computeJacobianhOver_q(args) { computeH(args); }   // old name, deprecated
 *  \endcode
 *
 *  - a class overriding the new name is called directly, through either name;
 *  - a class overriding the old name is called by Base::computeH;
 *  - with no override, Base::computeH -> Base::computeJacobianhOver_q -> Base::computeH,
 *    which then runs the default implementation (the flag is set).
 */
#ifndef RELATIONCOMPATIBILITY_HPP
#define RELATIONCOMPATIBILITY_HPP

/** Disable the deprecation warnings (GCC and clang), e.g. where a deprecated function is
 *  called on purpose to keep it working. */
#define SICONOS_IGNORE_DEPRECATED_BEGIN \
  _Pragma("GCC diagnostic push")        \
      _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#define SICONOS_IGNORE_DEPRECATED_END _Pragma("GCC diagnostic pop")

namespace siconos::modeling::detail {

/** RAII flag set while the deprecated name of a virtual function is being dispatched. */
class LegacyCallGuard {
 public:
  explicit LegacyCallGuard(bool& active) : active_(active) { active_ = true; }
  ~LegacyCallGuard() { active_ = false; }
  LegacyCallGuard(const LegacyCallGuard&) = delete;
  LegacyCallGuard& operator=(const LegacyCallGuard&) = delete;

 private:
  bool& active_;
};

}  // namespace siconos::modeling::detail

#endif  // RELATIONCOMPATIBILITY_HPP
