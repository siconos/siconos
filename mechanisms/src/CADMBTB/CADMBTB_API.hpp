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

/*! \addtogroup CADMBTB_API
 * \brief This module provides an API on the 3D modeler dedicated to the Multi
 * Bodies simulation.
 *
 *  It contains the API used by the module CADMBTB. <br>
 *  It provides the 3D modeler features used for the MBTB.
 *  @{
 */
#ifndef CADMBTBAPI
#define CADMBTBAPI
#include "TopoDS.hxx"

namespace siconos::mechanisms {

/** Updates the graphic.

    It assume that  CADMBTB_setGraphicContext has been called, else it does
    nothing.
*/
void CADMBTB_updateGraphic();

/** Initializes the CABMBTB library. It consists in allocating the working
 *  memory of n2qn1.
 *
 *  \param [in] NumberOfObj number of objects
 *  \param [in] NumberOfContacts number of contacts.
 */
void CADMBTB_init(unsigned int NumberOfObj, unsigned int NumberOfContacts);

/** CADMBTB_initContact
 *
 *  It is called to initialize the Bounding Box of the paremters. <br>
 *  It must be call as few as possible because of a bug in OCC: Computation
 *  becomming more and more slow. <br> It must be underline that the parameters
 *  used by qnb.f are not saved between the time steps.
 *
 *  \param contactId
 */
void CADMBTB_initContact(unsigned int contactId);

/** resets the  CABMBTB library: does nothing in current version */
void CADMBTB_reset();

/** Load a CAD file
 *
 * \param [in] id an identifier of the object (must be 0 =< id < NumberOfObj)
 * \param [in] fileName const char * , a CAD file
 */
void CADMBTB_loadCADFile(unsigned int id, const char* fileName);

/** Builds the graphical model of object, if it is not called, the object
 *  will be not draw in the 3D view.
 *
 *  \param [in] id identifier of the object (must be 0 =< id < NumberOfObj)
 */
void CADMBTB_buildGraphicalModel(unsigned int id);

/** Move an object from an other on. (useful for contact having the same
 *  position of DS) Implementation: It consists to apply the same DISPLACEMENT
 *  computed previously by the function CADMBTB_moveObjectFromQ. Warning :
 *  data::sStartTopoDS[id] is not updated
 *
 *  \param [in] idModel1 the identifier of the moved object
 *  \param [in] idModel2 the identifier of the referenced object
 */
void CADMBTB_moveModelFromModel(unsigned int idModel1, unsigned int idModel2);

/** Move a Graphical model from an object, indeed the graphical model has to
 *  follow the mechanical model. Implementation: It consists to set the current
 *  transformation using data::sGeomTrsf[idGraphicModel] computed in the function
 *  CADMBTB_moveObjectFromQ
 *
 *  \param [in] idGraphicModel the identifier of the graphical model
 *  \param [in] idModel the identifier of the referenced object
 */
void CADMBTB_moveGraphicalModelFromModel(unsigned int idGraphicModel,
                                         unsigned int idModel);

/** To move an object using quaternion.
 *  Implementation:
 *   1) It consists in computing the DISPLACEMENT beteween (x,y,z,q1,q2,q3,q4)
 * and the current position stored in data::sStartTopoDS[id]. 2) The displacement is
 * applyed 3) The data::sStartTopoDS[id] is updated
 *
 *   \param [in] id the identifier of the moved object
 *   \param [in] x translation in 1st direction
 *   \param [in] y translation in 2nd direction
 *   \param [in] z translation in 3rd direction
 *   \param[in] q1 first quaternion: cos(theta/2)
 *   \param [in] q2
 *   \param [in] q3
 *   \param [in] q4
 */
void CADMBTB_moveObjectFromQ(unsigned int id, double& x, double& y, double& z,
                             double& q1, double& q2, double& q3, double& q4);

/**  Set the location of an object WITHOUT moving it. Useful to defined the
 *   coordinate system attatched to an object during the initialization
 *
 *   \param[in] id the identifier of the object
 *   \param [in] x translation in 1st direction
 *   \param [in] y translation in 2nd direction
 *   \param [in] z translation in 3rd direction
 */
void CADMBTB_setLocation(unsigned int id, double& x, double& y, double& z);

/** Computes the bounding box (in parameter space), only once because OCC is
 * very slow
 *
 * \param id
 */
void CADMBTB_computeUVBounds(unsigned int id);

/** To get UV bound of the second elem (face or edge) of a shape
 *
 *  \param[in] id [in] Identifier of the shape.
 *  \param[out] U1  inf U bound.
 *  \param[out] U2  sup U bound.
 *  \param[out] V1  inf V bound.
 *  \param[out] V2  sup V bound.
 */
void CADMBTB_getUVBounds2(unsigned int id, double& U1, double& U2, double& V1,
                          double& V2);

/** To get UV bound of the first elem (face or edge) of a shape
 *
 *  \param[in]  id Identifier of the shape.
 *  \param[out] U1 inf U bound.
 *  \param[out] U2  sup U bound.
 *  \param[out] V1  inf V bound.
 *  \param[out] V2  sup V bound.
 */
void CADMBTB_getUVBounds(unsigned int id, double& U1, double& U2, double& V1,
                         double& V2);

/** To compute de distance between two objects, P1, P2 are the contact points in
 *  the abs frame. n is the nornmal, in the abs frame
 *
 *  \param[in] idContact id of contact(useful for drawing of artefacts)
 *  \param[in] id1 the identifier of the first object
 *  \param[in] id2 the identifier of the second object
 *  \param[out] X1 first coordinate of P1
 *  \param[out] Y1 second coordinate of P1
 *  \param[out] Z1 third coordinate of P1
 *  \param[out] X2 first coordinate of P2
 *  \param[out] Y2 second coordinate of P2
 *  \param[out] Z2 third coordinate of P2
 *  \param[out] nX first coordinate of n
 *  \param[out] nY second coordinate of n
 *  \param[out] nZ third coordinate of n
 *  \param normalFromFace1
 *  \param  MinDist distance
 */
void CADMBTB_getMinDistance(unsigned int idContact, unsigned int id1,
                            unsigned int id2, double& X1, double& Y1,
                            double& Z1, double& X2, double& Y2, double& Z2,
                            double& nX, double& nY, double& nZ,
                            unsigned int normalFromFace1, double& MinDist);

/** Declares the number of artefacts: ie graphical decoration(forces, normal,
 * P1P2)
 *
 *  \param nb number of artefacts
 */
void CADMBTB_setNbOfArtefacts(unsigned int nb);

/** To build a line artefact (P1P2)
 *
 * \param id
 * \param X1
 * \param Y1
 * \param Z1
 * \param X2
 * \param Y2
 * \param Z2
 */
void CADMBTB_buildLineArtefactLine(unsigned int id, double* X1, double* Y1,
                                   double* Z1, double* X2, double* Y2,
                                   double* Z2);

/** To build a oriented line artefact (n)
 *
 * \param id
 * \param X1
 * \param Y1
 * \param Z1
 * \param X2
 * \param Y2
 * \param Z2
 */
void CADMBTB_buildOrientedLineArtefactLine(unsigned int id, double* X1,
                                           double* Y1, double* Z1, double* X2,
                                           double* Y2, double* Z2);
/** To build a cylinder artefact (forces).
 *
 * \param id
 * \param X1
 * \param Y1
 * \param Z1
 * \param X2
 * \param Y2
 * \param Z2
 */
void CADMBTB_buildOrientedLineArtefactLine1(unsigned int id, double* X1,
                                            double* Y1, double* Z1, double* X2,
                                            double* Y2, double* Z2);

/** To build a cylinder artefact (forces).
 *
 * \param id
 * \param X1
 * \param Y1
 * \param Z1
 * \param X2
 * \param Y2
 * \param Z2
 * \param radius
 */
void CADMBTB_buildCylinderArtefactLine(unsigned int id, double* X1, double* Y1,
                                       double* Z1, double* X2, double* Y2,
                                       double* Z2, double* radius);

TopoDS_Shape CADMBTB_TopoDS(unsigned int num);

/*! @} */
}  // namespace siconos::mechanisms
#endif
