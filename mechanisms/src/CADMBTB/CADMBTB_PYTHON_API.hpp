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

/*! \addtogroup CADMBTB_PYTHON_API
  \brief This module contains the python API of the module CADMBTB.

  It provides the graphic functions allowing to plug an 3D view, from OCC. The
  3D view is optional.<br> If there is no 3D view, the internal graphical
  management will be ignored.<br> It also provides a functions to load an
  artefact CAD model and set internal parameters. <br> For more details see the
  documentation of the functions.
  *  @{
  */

#ifndef CADMBTBPYTHONAPI
#define CADMBTBPYTHONAPI

// OCC forwards
class AIS_InteractiveContext;
class V3d_View;

namespace siconos::mechanisms {

/** Load a step CAD file.
 *
 *  \param fileName a CAD file.
 *  \param trans transparency value
 */
void CADMBTB_loadArtefactCADFile(const char* fileName, double trans);

/** Initializes the graphical context
 *  (example : CADMBTB_setGraphicContext(display.Context) from PYTHON)
 *
 *  \param aisContext
 */
void CADMBTB_setGraphicContext(AIS_InteractiveContext& aisContext);

/** To set the current view.
 *
 *  \param aView
 */
void CADMBTB_setGraphicView(V3d_View& aView);

/** To disable graphic */
void CADMBTB_disableGraphic();

/** To set a double parameter.(extendable, without modifie the API)
 *  This type of function has been chosen to easely set any parameters without
 * modify the module API.
 *
 *  \param IdParam identifier of the param.<br> 0 for transparency.<br>
 *  \param idShape identifier of the shape.
 *  \param v : value.
 */
void CADMBTB_setShapeDParam(unsigned int IdParam, unsigned int idShape,
                            double v);

/** To set a double parameter.(extendable, without modifie the API)
 *  This type of function has been chosen to easily set any parameters without
 * modifying the module API.
 *
 *  \param IdParam identifier of the param.<br>
 *                        0 for transparency.<br>
 *  \param idContact identifier of the contact.
 *  \param idShape identifier of the shape of the contact (0 or 1).
 *  \param v value.
 */
void CADMBTB_setContactDParam(unsigned int IdParam, unsigned int idContact,
                              unsigned int idShape, double v);

/** To set a int parameter. (extendable, without modifie the API)
 *  This type of function has been chosen to easely set any parameters without
 * modifying the module API.
 *
 *  \param IdParam identifier of the param, used only for enable/disable the
 * dump of the 3D view. \param v value.
 */
void CADMBTB_setIParam(unsigned int IdParam, int v);

/** Manual dump of the 3D view. */
void CADMBTB_DumpGraphic();

/** Print of the distance
 *
 *  \param v
 */

void CADMBTB_print_dist(unsigned int v);

/*! @} */

}  // namespace siconos::mechanisms
#endif
