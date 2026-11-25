// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file Pylon.h
  Pylon class
*/
#pragma once

#include <nlohmann/json.hpp>

#include "SiconosVector.hpp"

namespace siconos::fem::cable {
/**   Description of a pylon (position, geometry ...)
 */
class Pylon {
 private:
  siconos::algebra::SiconosVector3 coordinates_ = {0., 0., 0.};

  /** True if it's a station pylon else false */
  bool isAStation_{false};

  /** radius (of the pulley for a station, of curvature in other cases) */
  double radius_{0.};

  /** distance between the pylon and the up-rope */
  double distanceToUpRope_{0.};

  /* distance between the pylon and the down-rope  */
  double distanceToDownRope_{0.};

  /** pylon height */
  double height_{0.};

  Pylon() = delete;

 public:
  Pylon(const Pylon&) = default;
  Pylon(Pylon&&) = default;
  Pylon& operator=(const Pylon&) = delete;
  Pylon& operator=(Pylon&&) = default;

  /** Build a pylon from a json input
    \param j json input
    \param is_station true if the pylon defines a station
    */
  Pylon(const nlohmann::json& j, bool is_station) : isAStation_(is_station) { from_json(j); };

  ~Pylon() noexcept = default;

  /** \return the coordinates of the pylon (read-only) */
  const siconos::algebra::SiconosVector3& coords() const { return coordinates_; }

  /** \return the  */
  inline auto get_radius() const { return radius_; }

  void from_json(const nlohmann::json& j);

  /** \return true if the pylon is a station (up or down) */
  bool isStation() const { return isAStation_; };

  /** Add the distance between up and down cables to the current pylon y value. */
  void shift_y();

  /** Screen display of the pylon parameters */
  void display() const;
};

// Required to be able to insert a Pylon into a set
// Comparison based on first coordinate (x)
bool operator<(const Pylon& p1, const Pylon& p2);

}  // namespace siconos::fem::cable
