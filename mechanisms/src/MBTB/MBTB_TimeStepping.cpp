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

#include "MBTB_TimeStepping.hpp"

#include "MBTB_PYTHON_API.hpp"    // For MBTB_updateDSFromSiconos
#include "MBTB_internalTool.hpp"  // For  _MBTB_updateContactFromDS();
// #define TS_DEBUG

void siconos::mechanisms::MBTB_TimeStepping::updateWorldFromDS() {
#ifdef TS_DEBUG
  printf("MBTB_TimeStepping::updateWordFromDS \n");
#endif
  MBTB_updateDSFromSiconos();
  mbtb::internal::MBTB_updateContactFromDS();
}
