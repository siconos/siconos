/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*!\file relay_cst.h
 * \brief Relay constants*/

#ifndef RELAY_CST_H
#define RELAY_CST_H

enum RELAY_SOLVER
{
  SICONOS_RELAY_PGS = 300,
  SICONOS_RELAY_ENUM = 301,
  SICONOS_RELAY_PATH = 302,
  SICONOS_RELAY_LEMKE = 303,
  SICONOS_RELAY_AVI_CAOFERRIS = 306,
  SICONOS_RELAY_AVI_CAOFERRIS_TEST = 307
};




extern const char* const   SICONOS_RELAY_PGS_STR;
extern const char* const   SICONOS_RELAY_PATH_STR;
extern const char* const   SICONOS_RELAY_ENUM_STR;
extern const char* const   SICONOS_RELAY_LEMKE_STR;
extern const char* const   SICONOS_RELAY_AVI_CAOFERRIS_STR;
extern const char* const   SICONOS_RELAY_AVI_CAOFERRIS_TEST_STR;
#endif
