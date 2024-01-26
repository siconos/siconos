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

/*! \file SiconosBulletShape.hpp
  \brief Siconos API to Bullet shapes
*/

#ifndef BULLETSICONOSSHAPE_H
#define BULLETSICONOSSHAPE_H

#include <BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>

#include <memory>
#include <vector>

// These classes are only used internaly in SiconosBulletCollisionManager.cpp
// They are not supposed to appear in user interface.

namespace siconos::collision::bullet::internal {

// We need a bit more space to hold mesh data
class SiconosMeshData : public btGImpactMeshShape {
 public:
  SiconosMeshData(btStridingMeshInterface *i)
      : btGImpactMeshShape(i), btScalarVertices(nullptr)
  {
  }
  ~SiconosMeshData() noexcept
  {
    if (btScalarVertices) delete[] btScalarVertices;
  }
  btScalar *btScalarVertices{nullptr};
  std::shared_ptr<btTriangleIndexVertexArray> btTriData{nullptr};
};

// Similarly, place to store height matrix
class SiconosHeightData : public btHeightfieldTerrainShape {
 public:
  SiconosHeightData(int width, std::shared_ptr<std::vector<btScalar>> data,
                    btScalar min_height, btScalar max_height)
      : btHeightfieldTerrainShape(width, data->size() / width, data->data(),
                                  0,  // scale ignored for PHY_FLOAT
                                  min_height, max_height, 2, PHY_FLOAT,
                                  false),  // up = z, flip = false
        _data(data),
        _min_height(min_height),
        _max_height(max_height)
  {
  }
  std::shared_ptr<std::vector<btScalar>> _data{nullptr};
  btScalar _min_height{0.}, _max_height{0.};
};

}  // namespace siconos::collision::bullet::internal
#endif
