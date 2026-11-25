// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SimulationTypeDef.hpp
 * \brief Typedef for simulation-related objects
 */

#ifndef SimulationTypes_H
#define SimulationTypes_H

namespace siconos::simulation {

enum class OsnspbType { DEFAULT, ED_SMOOTH_ACC, ED_IMPACT, ED_SMOOTH_POS };

enum SICONOS_OSNSP { SICONOS_OSNSP_DEFAULT = 0 };
enum SICONOS_OSNSP_ED {
  SICONOS_OSNSP_ED_SMOOTH_ACC,
  SICONOS_OSNSP_ED_IMPACT,
  SICONOS_OSNSP_ED_SMOOTH_POS
};
enum SICONOS_OSNSP_TS { SICONOS_OSNSP_TS_VELOCITY = 0, SICONOS_OSNSP_TS_POS = 1 };
constexpr int SICONOS_NB_OSNSP_TS = 1;
constexpr int SICONOS_NB_OSNSP_TSP = 2;

}  // namespace siconos::simulation
#endif
