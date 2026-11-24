// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SiconosBulletOptions.hpp
  \brief Simple classes to set Bullet options and get stats.
*/

#ifndef SICONOSBULLETOPTIONS_H
#define SICONOSBULLETOPTIONS_H

#include "SiconosSerialization.hpp"

namespace siconos::collision::bullet {

enum class SiconosBulletDimension { ThreeD = 0, TwoD = 1 };

/** Set of paramaters for Bullet API. All members have default values.*/
class SiconosBulletOptions {
 protected:
  ACCEPT_SERIALIZATION(SiconosBulletOptions);

 public:
  SiconosBulletDimension dimension{SiconosBulletDimension::ThreeD};
  double contactBreakingThreshold{0.02};
  double contactProcessingThreshold{0.03};
  double worldScale{1.0};
  bool useAxisSweep3{false};
  bool clearOverlappingPairCache{false};
  unsigned int perturbationIterations{3};
  unsigned int minimumPointsPerturbationThreshold{3};
  bool enableSatConvex{false};
  bool enablePolyhedralContactClipping{false};
  double Depth2D{0.04};
  double extrapolationCoefficient{0.};

  ~SiconosBulletOptions() noexcept = default;
};

/** To record stats during simu based on Bullet*/
class SiconosBulletStatistics {
 protected:
  ACCEPT_SERIALIZATION(SiconosBulletStatistics);

 public:
  int new_interactions_created{0};
  int existing_interactions_processed{0};
  int interaction_warnings{0};
  int interaction_destroyed{0};

  ~SiconosBulletStatistics() noexcept = default;
};

}  // namespace siconos::collision::bullet
#endif
