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
/*! \file TimeStepping.hpp
 *  \brief Time-Stepping simulation
 */

/*! \addtogroup CADMBTB_INTERNALTOOLS
   *  \brief This module contains the API of the internal funtions of the
   CADMBTB box.
   *
   *  It provides the geometrical functions. It consistes in computing the
   nearest points between objects. <br> It computes the cartesian coordinates
   and the corresponding normal. <br> For more details, see documentation about
   the functions.
   *  @{
   */
#ifndef CADMBTB_INTERNALTOOLS
#define CADMBTB_INTERNALTOOLS

#include <Standard_TypeDef.hxx>

namespace siconos::mechanisms {
//  /** To compute distance using OCC algorithm - Deprecated
//  *
//  *  \param[in] aFace1 the 1st object
//  *  \param[in] aFace2 the 2nd object
//  *  \param[out] X1 the 1st coordinate of the point attached to the 1st object
//  *  \param[out] Y1 the 2nd coordinate of the point attached to the 1st object
//  *  \param[out] Z1 the 3rd coordinate of the point attached to the 1st object
//  *  \param[out] X2 the 1st coordinate of the point attached to the 2nd object
//  *  \param[out] Y2 the 2nd coordinate of the point attached to the 2nd object
//  *  \param[out] Z2 the 3rd coordinate of the point attached to the 2nd object
//  *  \param[out] nX X component of the normal
//  *  \param[out] nY Y component of the normal
//  *  \param[out] nZ Z component of the normal
//  *  \param[out] MinDist the distance between the points
//  */
// void _CADMBTB_getMinDistanceFaceFace(const TopoDS_Face& aFace1, const
// TopoDS_Face& aFace2,
//                                      Standard_Real& X1, Standard_Real& Y1,
//                                      Standard_Real& Z1, Standard_Real& X2,
//                                      Standard_Real& Y2, Standard_Real& Z2,
//                                      Standard_Real& nX, Standard_Real& nY,
//                                      Standard_Real& nZ, Standard_Real&
//                                      MinDist);

// /** To compute distance using OCC algorithm - Deprecated
//  *
//  *  \param[in] idFace1 the 1st object
//  *  \param[in] idFace2 the 2nd object
//  *  \param[out] X1 the 1st coordinate of the point attached to the 1st object
//  *  \param[out] Y1 the 2nd coordinate of the point attached to the 1st object
//  *  \param[out] Z1 the 3rd coordinate of the point attached to the 1st object
//  *  \param[out] X2 the 1st coordinate of the point attached to the 2nd object
//  *  \param[out] Y2 the 2nd coordinate of the point attached to the 2nd object
//  *  \param[out] Z2 the 3rd coordinate of the point attached to the 2nd object
//  *  \param[out] nX  X component of the normal
//  *  \param[out] nY  Y component of the normal
//  *  \param[out] nZ  Z component of the normal
//  *  \param[out] MinDist the distance between the points
//  */
// void _CADMBTB_getMinDistanceFaceFace(unsigned int idFace1, unsigned int
// idFace2,
//                                      Standard_Real& X1, Standard_Real& Y1,
//                                      Standard_Real& Z1, Standard_Real& X2,
//                                      Standard_Real& Y2, Standard_Real& Z2,
//                                      Standard_Real& nX, Standard_Real& nY,
//                                      Standard_Real& nZ, Standard_Real&
//                                      MinDist);

/** To compute distance using n2qn1 algorithm This function manages the case
 * where the object idFace1 contains one or two faces
 *
 *  \param[in] idContact the identifier of the contact
 *  \param[in] idFace1 idFace1 the identifier of the 1st object
 *   containing one or two faces, redundancy parameter,
 *   must be equal to data::sNumberOfObj+(2*idContact-2*data::sNumberOfContacts)
 *  \param[in] idFace2 the identifier of the 2nd object containing one
 *   face, redundancy parameter, must be equal to
 *   data::sNumberOfObj+(2*idContact+1-2*data::sNumberOfContacts)
 *  \param[out] X1 the 1st coordinate of the point attached to the 1st object
 *  \param[out] Y1 the 2nd coordinate of the point attached to the 1st object
 *  \param[out] Z1 the 3rd coordinate of the point attached to the 1st object
 *  \param[out] X2 the 1st coordinate of the point attached to the 2nd object
 *  \param[out] Y2 the 2nd coordinate of the point attached to the 2nd object
 *  \param[out] Z2 the 3rd coordinate of the point attached to the 2nd object
 *  \param[in] normalFromFace1 if normalFromFace1 the normal is
 *   computed from the object idFace1
 *  \param[out] nX  X component of the normal
 *  \param[out] nY  Y component of the normal
 *  \param[out] nZ  Z component of the normal
 *  \param[out] MinDist the distance between the points
 */
void _CADMBTB_getMinDistanceFaceFace_using_n2qn1(
    unsigned int idContact, unsigned int idFace1, unsigned int idFace2,
    Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1, Standard_Real& X2,
    Standard_Real& Y2, Standard_Real& Z2, Standard_Real& nX, Standard_Real& nY,
    Standard_Real& nZ, unsigned int normalFromFace1, Standard_Real& MinDist);
/** To compute distance using n2qn1 algorithm
 * This function manages the case where the object idFace1 contains one or two
 faces
 *
 *  \param[in] idContact the identifier of the contact
 *  \param[in] idFace1 the identifier of the 1st object containing one or
 *   two faces, redundancy parameter, must be equal to
 *   data::sNumberOfObj+(2*idContact-2*data::sNumberOfContacts)
 *  \param[in] idFace2 the identifier of the 2nd object containing one
 *   edge, redundancy parameter, must be equal to
 *   data::sNumberOfObj+(2*idContact+1-2*data::sNumberOfContacts)
 *  \param[out] X1 the 1st coordinate of the point attached to the 1st object
 *  \param[out] Y1 the 2nd coordinate of the point attached to the 1st object
 *  \param[out] Z1 the 3rd coordinate of the point attached to the 1st object
 *  \param[out] X2 the 1st coordinate of the point attached to the 2nd object
 *  \param[out] Y2 the 2nd coordinate of the point attached to the 2nd object
 *  \param[out] Z2 the 3rd coordinate of the point attached to the 2nd object
 *  \param[out] nX  X component of the normal
 *  \param[out] nY  Y component of the normal
 *  \param[out] nZ  Z component of the normal
 *  \param normalFromFace1
 *  \param[out] MinDist the distance between the points
 */
void _CADMBTB_getMinDistanceFaceEdge_using_n2qn1(
    unsigned int idContact, unsigned int idFace1, unsigned int idFace2,
    Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1, Standard_Real& X2,
    Standard_Real& Y2, Standard_Real& Z2, Standard_Real& nX, Standard_Real& nY,
    Standard_Real& nZ, unsigned int normalFromFace1, Standard_Real& MinDist);
/*! @} */
}  // namespace siconos::mechanisms
#endif
