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

/**
 * \file BinaryCohesiveNSL.cpp
 * \brief Implementation of the BinaryCohesiveNSL cohesive zone model
 *
 * This file implements a binary (intact/broken) cohesive zone model where
 * the interface can be in one of two states:
 * - Intact (beta = 1): Can sustain traction up to sigma_c
 * - Broken (beta = 0): No cohesive traction, standard contact only
 *
 * \section sec_bczm_impl_damage Damage Evolution
 *
 * The damage parameter beta evolves based on the displacement jump delta:
 *
 * DOOR_SHAPE:
 * - If delta > delta_c: interface breaks instantaneously (beta = 0)
 * - If delta <= delta_c: interface remains intact (beta = 1)
 *
 * TRIANGLE_SHAPE:
 * - beta = min(beta_previous, 1 - delta/delta_c)
 * - Linear softening from beta=1 at delta=0 to beta=0 at delta=delta_c
 *
 * \section sec_bczm_impl_cohesion Cohesive Intensity Computation
 *
 * The cohesive intensity is computed as:
 * \f[ {cohesion} = -beta \cdot sigma_c \cdot surface \f] for the normal part
 * \f[ {cohesion} = -beta \cdot gamma \cdot sigma_c \cdot surface \cdot n \f] for the tangent
 * part
 *
 *
 * where:
 * - beta is the damage parameter
 * - sigma_c is the critical traction
 * - surface is the interface area
 * - gamma is the ration normal/tangent
 *
 * The negative sign indicates traction (pulling bodies together).
 *
 * \section sec_bczm_impl_internal Internal Variables Management
 *
 * The model stores 14 internal variables per interaction, organized as:
 * - BETA_SURFACE: [beta, surface_area]
 * - COHESION: 3D cohesive force vector
 * - DISPLACEMENT_JUMP: 3D displacement across interface
 * - COHESIVE_POINT_1/2: 3D contact points in global frame
 * - NORMAL/TANGENT_1/TANGENT_2: 3D local coordinate frame
 * - INITIAL_*: Initial values for persistency
 *
 * These variables are initialized in initializeInternalVariables() and
 * updated in updateInternalVariables() at each time step.
 *
 * \see BinaryCohesiveNSL.hpp for class documentation
 * \see CohesiveZoneModelNIFNSL for the base class
 */

#include "BinaryCohesiveNSL.hpp"

#include <algorithm>
#include <iostream>

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "NewtonEuler1DR.hpp"
#include "NewtonEulerR.hpp"
#include "RotationQuaternion.hpp"  // for orthoBaseFromVector
#include "SiconosVector.hpp"

// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT

#include "siconos_debug.h"

namespace siconos::mechanics::czm {

BinaryCohesiveNSL::BinaryCohesiveNSL(siconos::algebra::Index size)
    : siconos::modeling::CohesiveZoneModelNIFNSL(size) {}

BinaryCohesiveNSL::BinaryCohesiveNSL(double en, double et, double mu, double sigma_c,
                                     double delta_c, siconos::algebra::Index size)
    : siconos::modeling::CohesiveZoneModelNIFNSL(en, et, mu, size),
      _sigma_c(sigma_c),
      _delta_c(delta_c),
      _shape_type(ShapeType::DOOR_SHAPE) {}

BinaryCohesiveNSL::BinaryCohesiveNSL(double en, double et, double mu, double sigma_c,
                                     double delta_c, siconos::algebra::Index size,
                                     ShapeType shape_type)
    : siconos::modeling::CohesiveZoneModelNIFNSL(en, et, mu, size),
      _sigma_c(sigma_c),
      _delta_c(delta_c),
      _shape_type(shape_type) {
  if (_shape_type == ShapeType::TRIANGLE_SHAPE) {
    _slope = -1.0 / _delta_c;
  }
}

std::shared_ptr<siconos::algebra::blocks::SharedVector>
BinaryCohesiveNSL::initializeInternalVariables(siconos::modeling::Interaction& inter) {
  auto internalVariables_sp = std::make_shared<siconos::algebra::blocks::SharedVector>();
  internalVariables_sp->resize(BinaryCohesiveNSL::INTERNAL_VARIABLE_LENGTH);

  auto& internalVariables = *internalVariables_sp;

  auto rel = inter.relation();
  auto rel_NewtonEuler1DR = std::dynamic_pointer_cast<siconos::modeling::NewtonEuler1DR>(rel);

  if (rel_NewtonEuler1DR) {
    /* internalVariables(0) --> beta and surface */
    /* internalVariables(1) --> cohesion */
    /* internalVariables(2) --> displacement_jump */
    /* etc. */

    internalVariables[BinaryCohesiveNSL::COHESION] =
        std::make_shared<siconos::algebra::SiconosVector>(3);
    internalVariables[BinaryCohesiveNSL::DISPLACEMENT_JUMP] =
        std::make_shared<siconos::algebra::SiconosVector>(3);
    internalVariables[BinaryCohesiveNSL::BETA_SURFACE] =
        std::make_shared<siconos::algebra::SiconosVector>(2);

    // Initial value of beta = 1.0 (intact)
    (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])(0) = 1.0;
    // Initial value of surface = 1.0 (should be fixed correctly based on geometry)
    (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])(1) = 1.0;

    DEBUG_EXPR(std::cout << "\n The relation is of type NewtonEuler1DR" << std::endl;);

    // Store initial relative contact points
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_1] =
        std::make_shared<siconos::algebra::SiconosVector>(rel_NewtonEuler1DR->relPc1());
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_2] =
        std::make_shared<siconos::algebra::SiconosVector>(rel_NewtonEuler1DR->relPc2());

    // Store initial relative normal
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_NORMAL] =
        std::make_shared<siconos::algebra::SiconosVector>(rel_NewtonEuler1DR->relNc());

    // Compute tangent vectors from normal
    const auto& r_nc = rel_NewtonEuler1DR->relNc();
    siconos::algebra::SiconosVector3 t1, t2;

    if (!siconos::geometry::orthoBaseFromVector(
            const_cast<siconos::algebra::SiconosVector3&>(r_nc), t1, t2)) {
      THROW_EXCEPTION(
          "BinaryCohesiveNSL::initializeInternalVariables. Problem in calling "
          "orthoBaseFromVector");
    }

    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_1] =
        std::make_shared<siconos::algebra::SiconosVector>(t1);
    internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_2] =
        std::make_shared<siconos::algebra::SiconosVector>(t2);

    // Store current (absolute) contact points and normal
    internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_1] =
        std::make_shared<siconos::algebra::SiconosVector>(rel_NewtonEuler1DR->pc1());
    internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_2] =
        std::make_shared<siconos::algebra::SiconosVector>(rel_NewtonEuler1DR->pc2());
    internalVariables[BinaryCohesiveNSL::NORMAL] =
        std::make_shared<siconos::algebra::SiconosVector>(rel_NewtonEuler1DR->nc());

    // Compute tangent vectors from absolute normal
    const auto& nc = rel_NewtonEuler1DR->nc();
    siconos::algebra::SiconosVector3 abs_t1, abs_t2;

    if (!siconos::geometry::orthoBaseFromVector(
            const_cast<siconos::algebra::SiconosVector3&>(nc), abs_t1, abs_t2)) {
      THROW_EXCEPTION(
          "BinaryCohesiveNSL::initializeInternalVariables. Problem in calling "
          "orthoBaseFromVector for absolute frame");
    }

    internalVariables[BinaryCohesiveNSL::TANGENT_1] =
        std::make_shared<siconos::algebra::SiconosVector>(abs_t1);
    internalVariables[BinaryCohesiveNSL::TANGENT_2] =
        std::make_shared<siconos::algebra::SiconosVector>(abs_t2);

    // Compute initial displacement jump
    siconos::algebra::SiconosVector3 displacement_jump =
        rel_NewtonEuler1DR->pc2() - rel_NewtonEuler1DR->pc1();

    DEBUG_EXPR(siconos::algebra::print(displacement_jump););
    internalVariables[BinaryCohesiveNSL::INITIAL_DISPLACEMENT_JUMP] =
        std::make_shared<siconos::algebra::SiconosVector>(displacement_jump);

    DEBUG_EXPR(for (const auto& v : internalVariables) {
      if (v) siconos::algebra::print(*v);
    };);

  } else {
    // Simplified initialization for non-NewtonEuler1DR relations
    internalVariables[BinaryCohesiveNSL::COHESION] =
        std::make_shared<siconos::algebra::SiconosVector>(3);
    internalVariables[BinaryCohesiveNSL::BETA_SURFACE] =
        std::make_shared<siconos::algebra::SiconosVector>(2);

    (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])(0) = 1.0;  // beta
    (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])(1) = 1.0;  // surface
  }

  return internalVariables_sp;
}

void BinaryCohesiveNSL::updateInternalVariables(siconos::modeling::Interaction& inter) {
  DEBUG_BEGIN("void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)\n");

  auto internalVars = inter.internalVariables();
  auto internalVars_k = inter.internalVariables_k();

  if (!internalVars || !internalVars_k) {
    THROW_EXCEPTION(
        "BinaryCohesiveNSL::updateInternalVariables: internal variables not initialized");
  }

  auto& internalVariables = *internalVars;
  auto& internalVariables_k = *internalVars_k;

  double* beta = &(*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])(0);
  double* surface = &(*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])(1);
  double beta_k = (*internalVariables_k[BinaryCohesiveNSL::BETA_SURFACE])(0);

  double u_N = 0.0;
  double u_T = 0.0;
  double u_S = 0.0;

  if (beta_k > 0.0) {
    double delta = 0.0;

    auto rel = inter.relation();
    auto rel_NewtonEuler1DR =
        std::dynamic_pointer_cast<siconos::modeling::NewtonEuler1DR>(rel);

    if (rel_NewtonEuler1DR) {
      // Get stored initial configuration
      const auto& r_pc1_0 =
          *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_1];
      const auto& r_pc2_0 =
          *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_COHESIVE_POINT_2];
      const auto& r_nc_0 = *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_NORMAL];
      const auto& r_t1_0 = *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_1];
      const auto& r_t2_0 = *internalVariables[BinaryCohesiveNSL::INITIAL_RELATIVE_TANGENT_2];
      const auto& pos_0 = *internalVariables[BinaryCohesiveNSL::DISPLACEMENT_JUMP];

      // Access DS position from interaction
      const auto& ds_vars = inter.read_dynamical_systems_variables();
      // Note: In the modern API, we need to get q from DS differently
      // For now, use the relation's stored contact points

      auto& pc1_0 = *internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_1];
      auto& pc2_0 = *internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_2];
      auto& nc_0 = *internalVariables[BinaryCohesiveNSL::NORMAL];
      auto& t1_0 = *internalVariables[BinaryCohesiveNSL::TANGENT_1];
      auto& t2_0 = *internalVariables[BinaryCohesiveNSL::TANGENT_2];

      // Update contact points from current configuration
      // Note: In the old API this used computeContactPointsFromRelativeContactPoints
      // which needs to be adapted to the modern API. For now, we use the stored values.

      DEBUG_EXPR(std::cout << "pc1_0 is: " << pc1_0.transpose() << std::endl;
                 std::cout << "pc2_0 is: " << pc2_0.transpose() << std::endl;
                 std::cout << "nc_0 is: " << nc_0.transpose() << std::endl;);

      // Compute displacement
      siconos::algebra::SiconosVector3 pos = pc1_0 - pc2_0;

      DEBUG_EXPR(std::cout << "pos :" << pos.transpose() << std::endl;);
      siconos::algebra::SiconosVector3 u = pos - pos_0;

      u_N = u.dot(nc_0);
      u_T = u.dot(t1_0);
      u_S = u.dot(t2_0);

      DEBUG_EXPR(std::cout << "displacement jump u :" << u.transpose() << std::endl;);
      delta = u.norm();
      DEBUG_EXPR(std::cout << "delta :" << delta << std::endl;);

    } else {
      // Simplified case for non-NewtonEuler1DR
      auto size = this->size();
      if (size == 1) {
        delta = (*inter.y(0))(0);
      } else {
        delta = inter.y(0)->norm();
        const auto& y = *inter.y(0);
        u_T = y(1);
        u_S = y(2);
      }
    }

    /* Update beta (damage parameter) */
    DEBUG_PRINTF("beta = %e\n", *beta);
    DEBUG_PRINTF("beta_k = %e\n", beta_k);
    DEBUG_PRINTF("delta = %e\n", delta);

    if (_shape_type == ShapeType::DOOR_SHAPE) {
      if (delta > _delta_c) {
        DEBUG_PRINT("the interface is broken\n");
        *beta = 0.0;
      } else if (delta <= _delta_c && beta_k == 1.0) {
        DEBUG_PRINT("the interface is intact\n");
        *beta = 1.0;
      }
    } else if (_shape_type == ShapeType::TRIANGLE_SHAPE) {
      *beta = std::min(beta_k, 1.0 + _slope * delta);
      *beta = std::max(0.0, *beta);
    }
  } else {
    *beta = beta_k;
  }

  DEBUG_PRINTF("beta = %e\n", *beta);

  // Compute cohesion intensity
  double* cohesion = internalVariables[BinaryCohesiveNSL::COHESION]->data();

  // Initialize to zero
  for (int k = 1; k < this->size(); k++) {
    cohesion[k] = 0.0;
  }
  // Normal cohesion force (negative for traction)
  cohesion[0] = (*beta) * _sigma_c * (*surface);
  DEBUG_PRINTF("normal  cohesion intensity %4.2e\n", cohesion[0]);

  cohesion[1] = (*beta) * _gamma * _sigma_c * (*surface);
  DEBUG_PRINTF("tangent  cohesion intensity %4.2e\n", cohesion[1]);

  DEBUG_EXPR(siconos::algebra::print(*internalVariables[BinaryCohesiveNSL::COHESION]));

  DEBUG_END("void BinaryCohesiveNSL::updateInternalVariables(Interaction& inter)\n");
}

bool BinaryCohesiveNSL::isActiveAtLevel(siconos::modeling::Interaction& inter,
                                        unsigned int level) {
  auto internalVars = inter.internalVariables();
  if (!internalVars) return (level == 1);  // Default behavior

  double beta = (*(*internalVars)[BinaryCohesiveNSL::BETA_SURFACE])(0);

  // If interface is broken (beta == 0), use standard behavior
  // If intact, active at velocity level (level 1)
  if (beta > 0.0) {
    return (level == 1);
  }
  return (level == 1);
}

double* BinaryCohesiveNSL::cohesion(siconos::modeling::Interaction& inter) const {
  auto internalVars = inter.internalVariables();
  if (!internalVars) return nullptr;
  return (*internalVars)[BinaryCohesiveNSL::COHESION]->data();
}

double BinaryCohesiveNSL::beta(siconos::modeling::Interaction& inter) const {
  auto internalVars = inter.internalVariables();
  if (!internalVars) return 0.0;
  return (*(*internalVars)[BinaryCohesiveNSL::BETA_SURFACE])(0);
}

void BinaryCohesiveNSL::display() const {
  siconos::modeling::CohesiveZoneModelNIFNSL::display();
  std::cout << "=== BinaryCohesiveNSL data display ==============================="
            << std::endl;
  std::cout << "sigma_c: " << _sigma_c << std::endl;
  std::cout << "delta_c: " << _delta_c << std::endl;
  std::cout << "shape_type: " << (_shape_type == ShapeType::DOOR_SHAPE ? "DOOR" : "TRIANGLE")
            << std::endl;
  std::cout << "=================================================================="
            << std::endl;
}

void BinaryCohesiveNSL::displayInternalVariables(
    siconos::algebra::blocks::SharedVector& internalVariables) {
  std::cout << "=== BinaryCohesiveNSL Internal Variables ========================="
            << std::endl;

  // auto& internalVariables = *internalVars;

  // BETA_SURFACE
  if (internalVariables[BinaryCohesiveNSL::BETA_SURFACE]) {
    double beta = (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])(0);
    double surface = (*internalVariables[BinaryCohesiveNSL::BETA_SURFACE])(1);
    std::cout << "BETA (damage parameter): " << beta << " (1.0=intact, 0.0=broken)"
              << std::endl;
    std::cout << "Surface area: " << surface << std::endl;
  }

  // R_COHESION
  if (internalVariables[BinaryCohesiveNSL::COHESION]) {
    auto& coh = *internalVariables[BinaryCohesiveNSL::COHESION];
    std::cout << "Cohesion intensity COHESION: [" << coh(0);
    for (int i = 1; i < coh.size(); ++i) {
      std::cout << ", " << coh(i);
    }
    std::cout << "]" << std::endl;
  }

  // DISPLACEMENT_JUMP
  if (internalVariables[BinaryCohesiveNSL::DISPLACEMENT_JUMP]) {
    auto& disp = *internalVariables[BinaryCohesiveNSL::DISPLACEMENT_JUMP];
    std::cout << "Displacement jump: [" << disp(0);
    for (int i = 1; i < disp.size(); ++i) {
      std::cout << ", " << disp(i);
    }
    std::cout << "]" << std::endl;
  }

  // COHESIVE_POINT_1 and COHESIVE_POINT_2
  if (internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_1] &&
      internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_2]) {
    auto& pc1 = *internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_1];
    auto& pc2 = *internalVariables[BinaryCohesiveNSL::COHESIVE_POINT_2];
    std::cout << "Cohesive point 1: [" << pc1(0) << ", " << pc1(1) << ", " << pc1(2) << "]"
              << std::endl;
    std::cout << "Cohesive point 2: [" << pc2(0) << ", " << pc2(1) << ", " << pc2(2) << "]"
              << std::endl;
  }

  // NORMAL
  if (internalVariables[BinaryCohesiveNSL::NORMAL]) {
    auto& normal = *internalVariables[BinaryCohesiveNSL::NORMAL];
    std::cout << "Normal: [" << normal(0) << ", " << normal(1) << ", " << normal(2) << "]"
              << std::endl;
  }

  // TANGENT_1 and TANGENT_2
  if (internalVariables[BinaryCohesiveNSL::TANGENT_1] &&
      internalVariables[BinaryCohesiveNSL::TANGENT_2]) {
    auto& t1 = *internalVariables[BinaryCohesiveNSL::TANGENT_1];
    auto& t2 = *internalVariables[BinaryCohesiveNSL::TANGENT_2];
    std::cout << "Tangent 1: [" << t1(0) << ", " << t1(1) << ", " << t1(2) << "]" << std::endl;
    std::cout << "Tangent 2: [" << t2(0) << ", " << t2(1) << ", " << t2(2) << "]" << std::endl;
  }

  // INITIAL_DISPLACEMENT_JUMP
  if (internalVariables[BinaryCohesiveNSL::INITIAL_DISPLACEMENT_JUMP]) {
    auto& init_disp = *internalVariables[BinaryCohesiveNSL::INITIAL_DISPLACEMENT_JUMP];
    std::cout << "Initial displacement jump: [" << init_disp(0);
    for (int i = 1; i < init_disp.size(); ++i) {
      std::cout << ", " << init_disp(i);
    }
    std::cout << "]" << std::endl;
  }

  std::cout << "=================================================================="
            << std::endl;
}

}  // namespace siconos::mechanics::czm
