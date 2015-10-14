/*! \addtogroup CADMBTB_DATA
   * \brief This file contains the static variables used by the module CADMBTB. <br>
   *
   * About the memory management, the total number of objects is define by NB_OBJ. NB_OBJ is used to allocate the necessary memory. <br>
   * Each objet is typed with the CADMBTB_TYPE.
   *  @{
   */
#ifndef CADMBTBDATA
#define CADMBTBDATA


#include "Geom_Transformation.hxx"
#include "AIS_Shape.hxx"
#include "gp_Ax3.hxx"
#include "gp_Trsf.hxx"
#include "TopoDS_Face.hxx"
#include "TopExp_Explorer.hxx"
#include "V3d_View.hxx"

//! Maximal number of objects managed by the CADMBTB mudule.
#define NB_OBJ 200

//! The debug mode.
#define DEBUG_CONTACT
//! The type of object.
enum CADMBTB_TYPE
{
  CADMBTB_TYPE_NONE=0,
  CADMBTB_TYPE_FACE=1,
  CADMBTB_TYPE_EDGE=2
};


//!Number of obj, including both contacts and objects.
extern unsigned int sNumberOfObj;
//! pointer on the first Shape loaded
/*!
  FOR MBTB :
  sTopoDS[0 to MBTB::sNbOfBodies-1] contain the CAD model of the bodies.
  sTopoDS[MBTB::sNbOfBodies to end] contain the CAD model of the contacts.
*/
extern TopoDS_Shape sTopoDS[];
//!type of obj (Face or Contact)
extern unsigned int sTopoDSType[];
//! parameters borne inf
extern double sTopoDSBinf[];
//! parameters borne inf of the second shape
extern double sTopoDSBinf2[];
//!parameters borne sup
extern double sTopoDSBsup[];
//!parameters borne sup of the second shape
extern double sTopoDSBsup2[];
//! Graphical objects
extern AIS_Shape* spAISToposDS[];
//! Graphical objects transparency
extern double spAISTrans[];
//! It is the current position of the object. Named start because it is applied a deplacement.
extern gp_Ax3  sStartTopoDS[];
//! It is the deplacement registred in the function CADMBTB_moveObjectFromQ, and applied to the related objects (graphic model, contact model)
extern gp_Trsf  sTrsfTopoDS[];
//! The location of the object.
extern Geom_Transformation  sGeomTrsf[];

//! The graphical context, gave by the user.
extern AIS_InteractiveContext * pAIS_InteractiveContext;
//! The 3D view, gave by the user.
extern V3d_View * pV3d_View;
//!For automatic dump.
extern int sCmpDump;
//!For automatic dump.
extern unsigned int sDumpGraphic;
//!For manual dump.
extern unsigned int sCmpDumpMan;

//! Number of artefacts
extern unsigned int sNumberOfArtefacts;
//! Graphical model of the artefacts
extern AIS_Shape * sAISArtefacts[];
//! The minimal lenght allowing to draw an artefact.
extern double sMinLineLength;
//! Working memory for n2qn1
extern double * sWorkD;
//! Working memory for n2qn1
extern int * sWorkInt;
//! Number of contacts.
extern int sNumberOfContacts;
//!The verbose mode for cadprint_dist
extern unsigned int sCADPrintDist;
#endif
/*! @} */
