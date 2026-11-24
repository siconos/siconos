// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef TestContactUTILS_hpp
#define TestContactUTILS_hpp

#include <string>

#include "SiconosBulletOptions.hpp"

struct BounceParams {
  bool trace;
  bool dynamic;
  double size;
  double mass;
  double position;
  double timestep;
  double insideMargin;
  double outsideMargin;
  std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions> options{
      std::make_shared<siconos::collision::bullet::SiconosBulletOptions>()};

  void dump() const {
    printf("  trace:              %s\n", trace ? "on" : "off");
    printf("  dynamic:            %s\n", trace ? "on" : "off");
    printf("  size:               %.3g\n", size);
    printf("  mass:               %.3g\n", mass);
    printf("  position:           %.3g\n", position);
    printf("  insideMargin:       %.3g\n", insideMargin);
    printf("  outsideMargin:      %.3g\n", outsideMargin);
    printf("  breakingThreshold:  %.3g\n", options->contactBreakingThreshold);
    printf("  worldScale:         %.3g\n", options->worldScale);
  }
};

struct BounceResult {
  double bounce_error_sum;
  double bounce_error[6];
  int n_bounce_error;
  double final_position;
  double final_position_std;
  int num_interactions;
  int num_interaction_warnings;
  int max_simultaneous_contacts;
  double avg_simultaneous_contacts;
  double displacement_on_first_contact;
};

#endif
