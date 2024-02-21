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
#include "OccBody.hpp"

#include <TopoDS_Shape.hxx>
#include <boost/math/quaternion.hpp>

#include "OccContactEdge.hpp"
#include "OccContactFace.hpp"
#include "OccContactShape.hpp"
#include "OccUtils.hpp"
#include "SiconosVector.hpp"

siconos::mechanics::occ::OccBody::OccBody(std::shared_ptr<siconos::algebra::SiconosVector> position,
                               std::shared_ptr<siconos::algebra::SiconosVector> velocity,
                               double mass,
                               std::shared_ptr<siconos::algebra::SiconosMatrix> inertia)
    : NewtonEulerDS(position, velocity, mass, inertia),
      _contactShapes(std::make_shared<ContactShape_vector>()),
      _shapes(std::make_shared<TopoDS_Shape_vector>()) {}

void siconos::mechanics::occ::OccBody::addContactShape(
    OccContactShapeV& shape, std::shared_ptr<siconos::algebra::SiconosVector> pos,
    std::shared_ptr<siconos::algebra::SiconosVector> ori, unsigned int group) {
  OffSet offset = {0, 0, 0, 1, 0, 0, 0};

  if (pos) {
    offset[0] = (*pos)(0);
    offset[1] = (*pos)(1);
    offset[2] = (*pos)(2);
  }
  if (ori) {
    offset[3] = (*ori)(0);
    offset[4] = (*ori)(1);
    offset[5] = (*ori)(2);
    offset[6] = (*ori)(3);
  }

  _contactShapes->push_back(std::make_tuple(std::ref(shape), offset, group));

  updateContactShapes();

  auto computeUVBounds = [](auto& obj) { obj->computeUVBounds(); };
  std::visit(computeUVBounds, shape);
}

void siconos::mechanics::occ::OccBody::addShape(std::shared_ptr<TopoDS_Shape> shape,
                                     std::shared_ptr<siconos::algebra::SiconosVector> pos,
                                     std::shared_ptr<siconos::algebra::SiconosVector> ori) {
  OffSet offset = {0, 0, 0, 1, 0, 0, 0};
  if (pos) {
    offset[0] = (*pos)(0);
    offset[1] = (*pos)(1);
    offset[2] = (*pos)(2);
  }
  if (ori) {
    offset[3] = (*ori)(0);
    offset[4] = (*ori)(1);
    offset[5] = (*ori)(2);
    offset[6] = (*ori)(3);
  }

  _shapes->push_back(std::make_tuple(shape, offset));

  updateShapes();
}

void siconos::mechanics::occ::OccBody::updateContactShapes() {
  boost::math::quaternion<double> q((*_q)(3), (*_q)(4), (*_q)(5), (*_q)(6));

  for (auto& cs : *_contactShapes) {
    auto offset = std::get<1>(cs);

    boost::math::quaternion<double> pv =
        boost::math::quaternion<double>(0, offset[0], offset[1], offset[2]);

    boost::math::quaternion<double> rv = q * pv * boost::math::conj(q);

    boost::math::quaternion<double> r =
        q * boost::math::quaternion<double>(offset[3], offset[4], offset[5], offset[6]);

    siconos::algebra::SiconosVector fp = siconos::algebra::SiconosVector{7};
    fp(0) = (*_q)(0) + rv.R_component_2();
    fp(1) = (*_q)(1) + rv.R_component_3();
    fp(2) = (*_q)(2) + rv.R_component_4();
    fp(3) = r.R_component_1();
    fp(4) = r.R_component_2();
    fp(5) = r.R_component_3();
    fp(6) = r.R_component_4();

    auto getdata = [](auto& obj) { return obj->shape(); };
    auto shape_data = std::visit(getdata, std::get<0>(cs));

    siconos::mechanics::occ::occ_move(*shape_data, fp);
  }
}

void siconos::mechanics::occ::OccBody::updateShapes() {
  boost::math::quaternion<double> q((*_q)(3), (*_q)(4), (*_q)(5), (*_q)(6));

  for (auto& cs : *_shapes) {
    auto offset = std::get<1>(cs);

    boost::math::quaternion<double> pv =
        boost::math::quaternion<double>(0, offset[0], offset[1], offset[2]);

    boost::math::quaternion<double> rv = q * pv * boost::math::conj(q);

    boost::math::quaternion<double> r =
        q * boost::math::quaternion<double>(offset[3], offset[4], offset[5], offset[6]);

    siconos::algebra::SiconosVector fp = siconos::algebra::SiconosVector{7};
    fp(0) = (*_q)(0) + rv.R_component_2();
    fp(1) = (*_q)(1) + rv.R_component_3();
    fp(2) = (*_q)(2) + rv.R_component_4();
    fp(3) = r.R_component_1();
    fp(4) = r.R_component_2();
    fp(5) = r.R_component_3();
    fp(6) = r.R_component_4();

    occ_move(*(std::get<0>(cs)), fp);
  }
}

const siconos::mechanics::occ::OccContactShapeV& siconos::mechanics::occ::OccBody::contactShape(
    unsigned int id) const {
  return std::get<0>((*_contactShapes)[id]);
}

const TopoDS_Shape& siconos::mechanics::occ::OccBody::shape(unsigned int id) const {
  return *std::get<0>((*_shapes)[id]);
}
