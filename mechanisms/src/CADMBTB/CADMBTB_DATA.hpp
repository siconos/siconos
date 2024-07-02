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

/*! \addtogroup CADMBTB_DATA
 * \brief This file contains the static variables used by the module CADMBTB.
 * <br>
 *
 * About the memory management, the total number of objects is define by NB_OBJ.
 * NB_OBJ is used to allocate the necessary memory. <br> Each objet is typed
 * with the CADMBTB_TYPE.
 *  @{
 */
#ifndef CADMBTBDATA
#define CADMBTBDATA

//! Maximal number of objects managed by the CADMBTB mudule.
#define NB_OBJ 200

//! The debug mode.
#define DEBUG_CONTACT

// Forward decl for OCCT
class AIS_Shape;
class TopoDS_Shape;
class gp_Ax3;
class gp_Trsf;
class Geom_Transformation;
class V3d_View;

namespace siconos::mechanisms::data {
//! The type of object.
enum class CADMBTB_Type { None = 0, Face = 1, Edge = 2 };

//! Number of obj, including both contacts and objects.
extern unsigned int sNumberOfObj;
//! pointer on the first Shape loaded
/*!
  FOR MBTB :
  sTopoDS[0 to MBTB::sNbOfBodies-1] contain the CAD model of the bodies.
  sTopoDS[MBTB::sNbOfBodies to end] contain the CAD model of the contacts.
*/
extern TopoDS_Shape sTopoDS[];
//! type of obj (Face or Contact)
extern CADMBTB_Type sTopoDSType[];
//! parameters borne inf
extern double sTopoDSBinf[];
//! parameters borne inf of the second shape
extern double sTopoDSBinf2[];
//! parameters borne sup
extern double sTopoDSBsup[];
//! parameters borne sup of the second shape
extern double sTopoDSBsup2[];
//! Graphical objects
extern AIS_Shape* spAISToposDS[];
//! Graphical objects transparency
extern double spAISTrans[];
//! It is the current position of the object. Named start because it is
//! applied a deplacement.
extern gp_Ax3 sStartTopoDS[];
//! It is the deplacement registred in the function CADMBTB_moveObjectFromQ,
//! and applied to the related objects (graphic model, contact model)
extern gp_Trsf sTrsfTopoDS[];
//! The location of the object.
extern Geom_Transformation sGeomTrsf[];

//! The graphical context, gave by the user.
extern AIS_InteractiveContext* pAIS_InteractiveContext;
//! The 3D view, gave by the user.
extern V3d_View* pV3d_View;
//! For automatic dump.
extern int sCmpDump;
//! For automatic dump.
extern unsigned int sDumpGraphic;
//! For manual dump.
extern unsigned int sCmpDumpMan;

//! Number of artefacts
extern unsigned int sNumberOfArtefacts;
//! Graphical model of the artefacts
extern AIS_Shape* sAISArtefacts[];
//! The minimal lenght allowing to draw an artefact.
extern double sMinLineLength;
//! Working memory for n2qn1
extern double* sWorkD;
//! Working memory for n2qn1
extern int* sWorkInt;
//! Number of contacts.
extern int sNumberOfContacts;
//! The verbose mode for cadprint_dist
extern unsigned int sCADPrintDist;
}  // namespace siconos::mechanisms::data
#endif
/*! @} */
