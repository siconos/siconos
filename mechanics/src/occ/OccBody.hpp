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
/** \file OccBody.hpp
    \brief A Siconos Newton Euler dynamical system with
    associated contact shapes
 */

#ifndef OccBody_hpp
#define OccBody_hpp

#include <TopoDS_Shape.hxx>

#include "NewtonEulerDS.hpp"  // Base class
#include "OccContactShape.hpp"

namespace siconos::algebra {
class SiconosVector;
class SiconosMatrix;
}  // namespace siconos::algebra

// From OpenCASCADE
class TopoDS_Shape;

namespace siconos::mechanics::occ {

// class OccContactShape;

class OccBody : public siconos::modeling::NewtonEulerDS {
 private:
  using OffSet = std::array<double, 7>;
  using ContactShape_vector = std::vector<std::tuple<OccContactShapeV&, OffSet, int>>;
  using TopoDS_Shape_vector = std::vector<std::tuple<std::shared_ptr<TopoDS_Shape>, OffSet>>;

  std::shared_ptr<ContactShape_vector> _contactShapes{nullptr};
  std::shared_ptr<TopoDS_Shape_vector> _shapes{nullptr};

 public:
  // using NewtonEulerDS::NewtonEulerDS;

  virtual ~OccBody() noexcept = default;

  //! Constructor from a minimum set of data.
  //  \param position initial coordinates of this DynamicalSystem.
  //  \param velocity initial velocity of this DynamicalSystem.
  //  \param mass the mass.
  //  \param inertia the inertia matrix.
  //
  OccBody(std::shared_ptr<siconos::algebra::SiconosVector> position,
          std::shared_ptr<siconos::algebra::SiconosVector> velocity, double mass,
          std::shared_ptr<siconos::algebra::SiconosMatrix> inertia);

  /** Association of a contact shape.
   * \param shape the contact shape.
   * \param position relative position (x, y, z).
   * \param orientation relative orientation quaternion w, x, y, z
   * \param group contact group default 0
   */
  void addContactShape(OccContactShapeV& shape,
                       std::shared_ptr<siconos::algebra::SiconosVector> position =
                           std::shared_ptr<siconos::algebra::SiconosVector>(),
                       std::shared_ptr<siconos::algebra::SiconosVector> orientation =
                           std::shared_ptr<siconos::algebra::SiconosVector>(),
                       unsigned int group = 0);

  /** Association of a shape without contact.
   * \param shape the shape
   * \param position relative position (x, y, z).
   * \param orientation relative orientation quaternion w, x, y, z
   */
  void addShape(std::shared_ptr<TopoDS_Shape> shape,
                std::shared_ptr<siconos::algebra::SiconosVector> position =
                    std::shared_ptr<siconos::algebra::SiconosVector>(),
                std::shared_ptr<siconos::algebra::SiconosVector> orientation =
                    std::shared_ptr<siconos::algebra::SiconosVector>());

  /** Update positions and orientations of contact shapes.
   */
  void updateContactShapes();

  /** Update positions and orientations of shapes.
   */
  void updateShapes();

  /** Get an associated contact shape by its rank of association.
   *  \param id the number of the shape.
   */
  const OccContactShapeV& contactShape(unsigned int id) const;

  /** Get an associated shape by its rank of association.
   *  \param id the number of the shape.
   */
  const TopoDS_Shape& shape(unsigned int id) const;
};
}  // namespace siconos::mechanics::occ

#endif
