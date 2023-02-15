/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

/*! \file SimulationTypeDef.hpp
 * \brief Typedef for simulation-related objects
 */

#ifndef SimulationTypes_H
#define SimulationTypes_H


namespace siconos::simulation{

enum class OsnspbType
  {
    DEFAULT,
    ED_SMOOTH_ACC,
    ED_IMPACT,
    ED_SMOOTH_POS
  };

  
enum SICONOS_OSNSP { SICONOS_OSNSP_DEFAULT = 0 };
enum SICONOS_OSNSP_ED {
  SICONOS_OSNSP_ED_SMOOTH_ACC,
  SICONOS_OSNSP_ED_IMPACT,
  SICONOS_OSNSP_ED_SMOOTH_POS
};
enum SICONOS_OSNSP_TS { SICONOS_OSNSP_TS_VELOCITY = 0, SICONOS_OSNSP_TS_POS = 1 };
constexpr int SICONOS_NB_OSNSP_TS = 1;
constexpr int SICONOS_NB_OSNSP_TSP = 2;

}
#endif
