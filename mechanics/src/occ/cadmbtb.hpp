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
#ifndef cadmbtb_hpp
#define cadmbtb_hpp

#include <Standard_TypeDef.hxx>
// OpenCascade forward declarations
class gp_Pnt;
class gp_Dir;
class TopoDS_Face;
class TopoDS_Edge;

namespace siconos::mechanics::occ {

class OccContactFace;
class OccContactEdge;

namespace cadmbtb {

gp_Pnt FacePoint(const TopoDS_Face& face, Standard_Real u, Standard_Real v);
gp_Pnt EdgePoint(const TopoDS_Edge& edge, Standard_Real u);
gp_Dir FaceNormal(const TopoDS_Face& face, Standard_Real u, Standard_Real v);

void distanceFaceFace(const OccContactFace& csh1, const OccContactFace& csh2,
                      Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                      Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                      Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                      Standard_Real& MinDist);

void distanceFaceEdge(const OccContactFace& sh1, const OccContactEdge& sh2, Standard_Real& X1,
                      Standard_Real& Y1, Standard_Real& Z1, Standard_Real& X2,
                      Standard_Real& Y2, Standard_Real& Z2, Standard_Real& nX,
                      Standard_Real& nY, Standard_Real& nZ, Standard_Real& MinDist);

}  // namespace cadmbtb
}  // namespace siconos::mechanics::occ
#endif
