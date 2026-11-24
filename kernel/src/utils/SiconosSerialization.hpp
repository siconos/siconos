/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*! \file SiconosSerialization.hpp
  serialization for Siconos
*/

#ifndef SiconosSerialization_hpp
#define SiconosSerialization_hpp

// tells include-what-you-use to keep this file
// and not to suggest boost or alike.
// IWYU pragma: begin_exports

namespace boost {
namespace serialization {
class access;
}
}  // namespace boost

/** install serialization hooks. Must be used inside a protected zone
    of class definition
    \parameter a class name
*/
#define ACCEPT_SERIALIZATION(CLASS)                             \
  typedef void serializable;                                    \
  template <typename Archive>                                   \
  friend void siconos_io(Archive&, CLASS&, const unsigned int); \
  friend class boost::serialization::access

// IWYU pragma: end_exports

#endif
