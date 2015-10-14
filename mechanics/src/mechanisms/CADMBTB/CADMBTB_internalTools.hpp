/*! \addtogroup CADMBTB_INTERNALTOOLS
   *  \brief This module contains the API of the internal funtions of the CADMBTB box.
   *
   *  It provides the geometrical functions. It consistes in computing the nearest points between objects. <br>
   It computes the cartesian coordinates and the corresponding normal. <br>
   For more details, see documentation about the functions.
   *  @{
   */
#ifndef CADMBTB_INTERNALTOOLS
#define CADMBTB_INTERNALTOOLS
#include "TopoDS_Face.hxx"
/** To compute distance using OCC algorithm. Deprecated.
 * \param  [in] aFace1 TopoDS_Face&  the first object.
 * \param  [in] aFace2 TopoDS_Face&  the second object.
 * \param  [out] X1 Standard_Real&  the first coordinate of the point attached to the first object.
 * \param  [out] Y1 Standard_Real&  the second coordinate of the point attached to the first object.
 * \param  [out] Z1 Standard_Real&  the third coordinate of the point attached to the first object.
 * \param  [out] X2 Standard_Real&  the first coordinate of the point attached to the second object.
 * \param  [out] Y2 Standard_Real&  the second coordinate of the point attached to the second object.
 * \param  [out] Z2 Standard_Real&  the third coordinate of the point attached to the second object.
 * \param  [out] nX  X component of the normal
 * \param  [out] nY  Y component of the normal
 * \param  [out] nZ  Z component of the normal
 * \param  [out] MinDist Standard_Real&  the distance between the points.
 */
void _CADMBTB_getMinDistanceFaceFace(const TopoDS_Face& aFace1, const TopoDS_Face& aFace2,
                                     Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                     Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                     Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                                     Standard_Real& MinDist);

/** To compute distance using OCC algorithm. Deprecated.
 * \param  [in] idFace1 TopoDS_Face&  the first object.
 * \param  [in] idFace2 TopoDS_Face&  the second object.
 * \param  [out] X1 Standard_Real&  the first coordinate of the point attached to the first object.
 * \param  [out] Y1 Standard_Real&  the second coordinate of the point attached to the first object.
 * \param  [out] Z1 Standard_Real&  the third coordinate of the point attached to the first object.
 * \param  [out] X2 Standard_Real&  the first coordinate of the point attached to the second object.
 * \param  [out] Y2 Standard_Real&  the second coordinate of the point attached to the second object.
 * \param  [out] Z2 Standard_Real&  the third coordinate of the point attached to the second object.
 * \param  [out] nX  X component of the normal
 * \param  [out] nY  Y component of the normal
 * \param  [out] nZ  Z component of the normal
 * \param  [out] MinDist Standard_Real&  the distance between the points.
 */
void _CADMBTB_getMinDistanceFaceFace(unsigned int idFace1, unsigned int idFace2,
                                     Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                     Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                     Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                                     Standard_Real& MinDist);
/** To compute distance using n2qn1 algorithm. This function manages the case where the object idFace1 contains one or two faces.
 * \param  [in] idContact unsigned int  the identifier of the contact.
 * \param  [in] idFace1 unsigned int idFace1 the identifier of the first object containing one or two faces, redundancy parameter, must be equal to sNumberOfObj+(2*idContact-2*sNumberOfContacts).
 * \param  [in] idFace2 unsigned int  the identifier of the second object containing one face, redundancy parameter, must be equal to sNumberOfObj+(2*idContact+1-2*sNumberOfContacts).
 * \param  [out] X1 Standard_Real&  the first coordinate of the point attached to the first object.
 * \param  [out] Y1 Standard_Real&  the second coordinate of the point attached to the first object.
 * \param  [out] Z1 Standard_Real&  the third coordinate of the point attached to the first object.
 * \param  [out] X2 Standard_Real&  the first coordinate of the point attached to the second object.
 * \param  [out] Y2 Standard_Real&  the second coordinate of the point attached to the second object.
 * \param  [out] Z2 Standard_Real&  the third coordinate of the point attached to the second object.
 * \param  [in] normalFromFace1 unsigned int , if normalFromFace1 the normal is computed from the object idFace1.
 * \param  [out] nX  X component of the normal
 * \param  [out] nY  Y component of the normal
 * \param  [out] nZ  Z component of the normal
 * \param  [out] MinDist Standard_Real&  the distance between the points.
 */
void _CADMBTB_getMinDistanceFaceFace_using_n2qn1(unsigned int idContact,unsigned int idFace1, unsigned int idFace2,
                                                 Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                                 Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                                 Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                                                 unsigned int normalFromFace1,
                                                 Standard_Real& MinDist);
/** To compute distance using n2qn1 algorithm.
 * This function manages the case where the object idFace1 contains one or two faces.
 * \param  [in] idContact unsigned int  the identifier of the contact.
 * \param  [in] idFace1 unsigned int  the identifier of the first object containing one or two faces, redundancy parameter, must be equal to sNumberOfObj+(2*idContact-2*sNumberOfContacts).
 * \param  [in] idFace2 unsigned int  the identifier of the second object containing one edge, redundancy parameter, must be equal to sNumberOfObj+(2*idContact+1-2*sNumberOfContacts).
 * \param  [out] X1 Standard_Real&  the first coordinate of the point attached to the first object.
 * \param  [out] Y1 Standard_Real&  the second coordinate of the point attached to the first object.
 * \param  [out] Z1 Standard_Real&  the third coordinate of the point attached to the first object.
 * \param  [out] X2 Standard_Real&  the first coordinate of the point attached to the second object.
 * \param  [out] Y2 Standard_Real&  the second coordinate of the point attached to the second object.
 * \param  [out] Z2 Standard_Real&  the third coordinate of the point attached to the second object.
 * \param  [out] nX  X component of the normal
 * \param  [out] nY  Y component of the normal
 * \param  [out] nZ  Z component of the normal
 * \param normalFromFace1
 * \param  [out] MinDist 
Standard_Real&  the distance between the points.
 */
void _CADMBTB_getMinDistanceFaceEdge_using_n2qn1(
  unsigned int idContact,unsigned int idFace1, unsigned int idFace2,
  Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
  Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
  Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
  unsigned int normalFromFace1, 
  Standard_Real& MinDist);
/*! @} */
#endif
