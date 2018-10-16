
/*! \addtogroup CADMBTB_API
 * \brief This module provides an API on the 3D modeler dedicated to the Multi Bodies simulation.
 *
 *  It contains the API used by the module CADMBTB. <br>
 *  It provides the 3D modeler features used for the MBTB.
 *  @{
 */
#ifndef CADMBTBAPI
#define CADMBTBAPI
#include "TopoDS.hxx"


//! It updates the graphic.
/*!
  It assume that  CADMBTB_setGraphicContext has been called, else it does nothing.
*/
void CADMBTB_updateGraphic();
//! It initializes the CABMBTB library. It consists in allocating the working memory of n2qn1.
/*!
  \param [in] NumberOfObj unsigned int .
  \param [in] NumberOfContacts number of contacts.
 */
void CADMBTB_init(unsigned int NumberOfObj,unsigned int NumberOfContacts);

/** CADMBTB_initContact
 * It is called to initialize the Bounding Box of the paremters. <br>
 * It must be call as few as possible because of a bug in OCC: Computation becomming more and more slow. <br>
 * It must be underline that the parameters used by qnb.f are not saved between the time steps.
 * \param contactId
 */
void CADMBTB_initContact(unsigned int contactId);
//! CADMBTB_reset
/*!
  It resets the  CABMBTB library: do nothing in current version.
 */
void CADMBTB_reset();

/** CADMBTB_loadCADFile.
 * Load a CAD file
 *  \param [in] id unsigned int  an identifier of the object (must be 0 =< id < NumberOfObj)
 *  \param [in] fileName const char * , a CAD file.
 */
void CADMBTB_loadCADFile(unsigned int id, const char * fileName);

/** CADMBTB_buildGraphicalModel
 * It builds the graphical model of object, if it is not called, the object will be not draw in the 3D view.
 * \param [in] id unsigned int  identifier of the object (must be 0 =< id < NumberOfObj)
 */
void CADMBTB_buildGraphicalModel(unsigned int id);

/** CADMBTB_moveModelFromModel
 * To move a object from an other on. (useful for contact having the same position of DS)
 *  Implementation:
 *  It consists to apply the same DISPLACEMENT computed previously by the function CADMBTB_moveObjectFromQ.
 *  Warning : sStartTopoDS[id] is not updated.
 *  \param [in] idModel1 unsigned int , the identifier of the moved object
 *  \param [in] idModel2 unsigned int idModel2, the identifier of the referenced object
 */
void CADMBTB_moveModelFromModel(unsigned int idModel1, unsigned int idModel2);

/** CADMBTB_moveGraphicalModelFromModel
 * To move a Graphical model from an object, indeed the graphical model has to follow the mechanical model.
 * Implementation:
 * It consists to set the current transformation using sGeomTrsf[idGraphicModel] computed in the function CADMBTB_moveObjectFromQ.
 * \param [in] idGraphicModel unsigned int , the identifier of the graphical model
 * \param [in] idModel unsigned int , the identifier of the referenced object
 */
void CADMBTB_moveGraphicalModelFromModel(unsigned int idGraphicModel, unsigned int idModel);

/** CADMBTB_moveObjectFromQ
  To move an object using quaternion.
  Implementation:
  1) It consists in computing the DISPLACEMENT beteween (x,y,z,q1,q2,q3,q4) and the current position stored in sStartTopoDS[id].
  2) The displacement is applayed
  3) The sStartTopoDS[id] is updated
  \param [in] id unsigned int , the identifier of the moved object
  \param [in] x double& x, x translation
  \param [in] y double& y, y translation
  \param [in] z double& z, z translation
  \param [in] q1 double& q1, first quaternion: cos(theta/2)
  \param [in] q2 double& q2,
  \param [in] q3 double& q3,
  \param [in] q4 double& q4,
 */
void CADMBTB_moveObjectFromQ(unsigned int id,double& x,double& y, double& z, double& q1,double& q2,double& q3,double& q4);

/** CADMBTB_setLocation
 * Set the location of an object WITHOUT moving it. Useful to defined the coordinate system attatched to an object during the initialization.
 * \param [in] id unsigned int idModel, the identifier of the object
 * \param [in] x double& x, x translation
 * \param [in] y double& y, y translation
 * \param [in] z double& z, z translation
 */
void CADMBTB_setLocation(unsigned int id, double& x,double& y, double& z);

/** CADMBTB_computeUVBounds
 *  It comuptes the boubing box (in parameter space), do only once because OCC is very slow.
 * \param id
 */
void CADMBTB_computeUVBounds(unsigned int id);

/** CADMBTB_getUVBounds2
 * To get UV bound of the second elem (face or edge) of a shape
 * \param id [in]: Identifier of the shape.
 * \param U1 [out]: inf U bound.
 * \param U2 [out]: sup U bound.
 * \param V1 [out]: inf V bound.
 * \param V2 [out]: sup V bound.
 */
void CADMBTB_getUVBounds2(unsigned int id, double& U1, double& U2, double& V1, double& V2);

/** CADMBTB_getUVBounds
 * To get UV bound of the first elem (face or edge) of a shape
 * \param id [in]: Identifier of the shape.
 * \param U1 [out]: inf U bound.
 * \param U2 [out]: sup U bound.
 * \param V1 [out]: inf V bound.
 * \param V2 [out]: sup V bound.
 */
void CADMBTB_getUVBounds(unsigned int id, double& U1, double& U2, double& V1, double& V2);

/** CADMBTB_getMinDistance
 * To compute de distance between two objects, P1, P2 are the contact points in the abs frame.
 * n is the nornmal, in the abs frame.
 * \param[in] idContact unsigned int : id of contact(useful for drawing of artefacts).
 * \param[in] id1 unsigned int id1: the identifier of the first object.
 * \param[in] id2 unsigned int id2: the identifier of the second object.
 * \param [out] X1 double first coordinate of P1.
 * \param [out] Y1 double second coordinate of P1.
 * \param [out] Z1 double third coordinate of P1.
 * \param [out] X2 double first coordinate of P2.
 * \param [out] Y2 double second coordinate of P2.
 * \param [out] Z2 double third coordinate of P2.
 * \param [out] nX double first coordinate of n.
 * \param [out] nY double second coordinate of n.
 * \param [out] nZ double third coordinate of n.
 * \param normalFromFace1
 * \param  MinDist double distance.
 */
void CADMBTB_getMinDistance
(unsigned int idContact,unsigned int id1, unsigned int id2,
 double& X1, double& Y1, double& Z1,
 double& X2, double& Y2, double& Z2,
 double& nX, double& nY, double& nZ, unsigned int normalFromFace1,
 double& MinDist);

/** CADMBTB_setNbOfArtefacts
 * To declare the number of artefacts: ie graphical decoration( forces, normal, P1P2)
 * \param nb number of artefacts
 */
void CADMBTB_setNbOfArtefacts(unsigned int nb);

/** CADMBTB_buildLineArtefactLine
 * To build a line artefact (P1P2).
 * \param id
 * \param X1
 * \param Y1
 * \param Z1
 * \param X2
 * \param Y2
 * \param Z2
 */
void CADMBTB_buildLineArtefactLine(unsigned int id,  double* X1, double* Y1, double* Z1,
                                   double* X2, double* Y2, double* Z2);

/** CADMBTB_buildOrientedLineArtefactLine
 *  To build a oriented line artefact (n).
 * \param id
 * \param X1
 * \param Y1
 * \param Z1
 * \param X2
 * \param Y2
 * \param Z2
 */
void CADMBTB_buildOrientedLineArtefactLine(unsigned int id,  double* X1, double* Y1, double* Z1,
    double* X2, double* Y2, double* Z2);
/** CADMBTB_buildCylinderArtefactLine
 * To build a cylinder artefact (forces).
 * \param id
 * \param X1
 * \param Y1
 * \param Z1
 * \param X2
 * \param Y2
 * \param Z2
 */
void CADMBTB_buildOrientedLineArtefactLine1(unsigned int id,  double* X1, double* Y1, double* Z1, double* X2, double* Y2, double* Z2);


/** CADMBTB_buildCylinderArtefactLine
 * To build a cylinder artefact (forces).
 * \param id
 * \param X1
 * \param Y1
 * \param Z1
 * \param X2
 * \param Y2
 * \param Z2
 * \param radius
 */
void CADMBTB_buildCylinderArtefactLine(unsigned int id,  double* X1, double* Y1, double* Z1,
                                       double* X2, double* Y2, double* Z2, double *radius);


TopoDS_Shape CADMBTB_TopoDS(unsigned int num);

/*! @} */

#endif
