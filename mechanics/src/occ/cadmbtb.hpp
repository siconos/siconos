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
#include <array>
#include <memory>
#include <span>

#include "ContactShapeDistance.hpp"

// OpenCascade forward declarations
class gp_Pnt;
class gp_Dir;
class TopoDS_Face;
class TopoDS_Edge;

namespace siconos::mechanics::occ {

class OccContactFace;
class OccContactEdge;

namespace cadmbtb {

// Functions required in mechanics APIS (User, Python ...)
auto distanceFaceFace(std::shared_ptr<OccContactFace> csh1,
                      std::shared_ptr<OccContactFace> csh2)
    -> ContactShapeDistance;

auto distanceFaceEdge(std::shared_ptr<OccContactFace> sh1,
                      std::shared_ptr<OccContactEdge> sh2)
    -> ContactShapeDistance;

namespace tools {

// template <std::size_t N>
// auto wrap_n2qn1(const std::array<double, N>& binf, const
// std::array<double, N>& bsup,
//                 const TopoDS_Face& face1, const TopoDS_Edge& edge2);

// Tools, functions, with direct calls to OCC functions
// - Supposed to be used only internally (mechanics or mechanisms)
// - Depends only on OpenCASCADE objects
// - Works with OCC 7.7

/** computes distance between two faces and its gradient

    \param x array with parameters on both surfaces (u1,v1,u2,v2)
    \param face1 occt first face
    \param face2 occt second face
    \return a tuple, first elem = square magnitude of distance func,
   second = array with gradient components
 */
auto myf_FaceFace(std::span<const double> x, std::shared_ptr<TopoDS_Face> face1,
                  std::shared_ptr<TopoDS_Face> face2)
    -> std::tuple<double, std::unique_ptr<std::array<double, 4>>>;

/** Computes distance between two faces and its gradient

      \param[in] x array with parameters on surface and edge [u_s, v_s,
   u_e] \param[in] face occt face \param[in] edge occt edge \return  a
   tuple, first elem = square magnitude of distance func, second = array
   with gradient components
*/
auto myf_FaceEdge(std::span<const double> x, std::shared_ptr<TopoDS_Face> face,
                  std::shared_ptr<TopoDS_Edge> edge)
    -> std::tuple<double, std::unique_ptr<std::array<double, 3>>>;

/** \return a point from the parameters values on the face

      \param[in] face occ face
      \param[in] u first parameter
      \param[in] v second parameter
*/
gp_Pnt FacePoint(const TopoDS_Face& face, Standard_Real u, Standard_Real v);

/** \return  a point from the parameter value on the Edge

      \param[in] edge occ edge
      \param[in] u
*/
gp_Pnt EdgePoint(const TopoDS_Edge& edge, Standard_Real u);

/** \return the normal at a point given from the parameters on the face.

      \param[in] face occ face
      \param[in] u first parameter
      \param[in] v second parameter
*/
gp_Dir FaceNormal(const TopoDS_Face& face, Standard_Real u, Standard_Real v);

}  // namespace tools

}  // namespace cadmbtb
}  // namespace siconos::mechanics::occ
#endif
