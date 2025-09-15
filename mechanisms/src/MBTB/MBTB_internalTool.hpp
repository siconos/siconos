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

/*! \addtogroup MBTB_INTERNAL_TOOL
   \brief This file contains the internal tools of the MBTB.

   * It consists in updating the CADMBTB from the simulation step. <br>
   * It also manages the output data.
   *  @{
   */
#ifndef INTERNALTOOLMBTB
#define INTERNALTOOLMBTB
#include <stdio.h>

#include <string>
#define PRINT_FORCE_CONTACTS
#define MBTB_PRINT_DIST
//! It updates the contacts CAD model from the body.
/*!

 */
void _MBTB_updateContactFromDS();

/** It updates the contacts CAD model from the body.
 * \param [in] numDS int,  update the cad model of contact related to the ds of id numDS.
 */
void _MBTB_updateContactFromDS(int numDS);

FILE *_MBTB_open(std::string filename, std::string args);

void _MBTB_close(FILE *);

/**!It prints the header of the output file.
 * \param fp output file
 */
void _MBTB_printHeader(FILE *fp);

/**It prints the current state in the output file.
 * \param fp output file
 */
void _MBTB_printStep(FILE *fp);

/** It displays the current state on std output.
 */
void _MBTB_displayStep();

/** It performs a step including the siconos call and the graphical update.
 */
void _MBTB_STEP();

#endif
/*! @} */
