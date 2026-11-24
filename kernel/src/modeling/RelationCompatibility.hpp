// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file RelationCompatibility.hpp
 *  \brief Helpers to keep the deprecated names of the relation API working.
 *
 *  When a virtual function of a relation is renamed (e.g. LagrangianScleronomousR::
 *  computeJacobianhOver_q -> computeContactMatrix), a user class that still overrides the old
 * name must keep being called. The pattern is:
 *
 *  \code
 *  void Base::computeContactMatrix(args) {          // new name, the one called by the kernel
 *    if (!computeContactMatrixLegacy_) {            // first entry: give a chance to an old
 * override LegacyCallGuard guard(computeContactMatrixLegacy_); SICONOS_IGNORE_DEPRECATED_BEGIN
 *      computeJacobianhOver_q(args);    // virtual: user override, or Base version below
 *      SICONOS_IGNORE_DEPRECATED_END
 *      return;
 *    }
 *    ...default implementation...
 *  }
 *  void Base::computeJacobianhOver_q(args) { computeContactMatrix(args); }   // old name,
 * deprecated
 *  \endcode
 *
 *  - a class overriding the new name is called directly, through either name;
 *  - a class overriding the old name is called by Base::computeContactMatrix;
 *  - with no override, Base::computeContactMatrix -> Base::computeJacobianhOver_q ->
 * Base::computeContactMatrix, which then runs the default implementation (the flag is set).
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
