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

#include "SiconosShape.hpp"

#include "ShapeVisitor.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::collision::SiconosBox::SiconosBox(double width, double height, double depth)
    : SiconosShape(), _dimensions(std::make_shared<siconos::algebra::SiconosVector>(3)) {
  (*_dimensions)(0) = width;
  (*_dimensions)(1) = height;
  (*_dimensions)(2) = depth;
}

void siconos::collision::SiconosBox::setDimensions(double width, double height, double depth) {
  (*_dimensions)(0) = width;
  (*_dimensions)(1) = height;
  (*_dimensions)(2) = depth;
  _version++;
}

void siconos::collision::SiconosBox::setDimensions(
    const siconos::algebra::SiconosVector& dim) {
  (*_dimensions)(0) = dim(0);
  (*_dimensions)(1) = dim(1);
  (*_dimensions)(2) = dim(2);
  _version++;
}

siconos::collision::SiconosConvexHull::SiconosConvexHull(
    std::shared_ptr<siconos::algebra::SiconosMatrix> vertices)
    : SiconosShape(), _vertices(vertices) {
  if (_vertices && _vertices->cols() != 3)
    THROW_EXCEPTION("Convex hull vertices matrix must have 3 columns.");
}

siconos::collision::SiconosConvexHull::SiconosConvexHull(
    Eigen::Ref<siconos::algebra::SiconosMatrix> vertices)
    : SiconosShape(), _vertices(std::make_shared<siconos::algebra::SiconosMatrix>(vertices)) {
  if (_vertices && _vertices->cols() != 3)
    THROW_EXCEPTION("Convex hull vertices matrix must have 3 columns.");
}

siconos::collision::SiconosMesh::SiconosMesh(
    std::shared_ptr<std::vector<unsigned int>> indexes,
    std::shared_ptr<siconos::algebra::SiconosMatrix> vertices)
    : SiconosShape(), _indexes(indexes), _vertices(vertices) {
  if (!_indexes || (_indexes->size() % 3) != 0)
    THROW_EXCEPTION("Mesh indexes size must be divisible by 3.");
  if (!_vertices || _vertices->rows() != 3)
    THROW_EXCEPTION("Mesh vertices matrix must have 3 columns.");
}
siconos::collision::SiconosMesh::SiconosMesh(
    const std::vector<unsigned int>& indexes,
    Eigen::Ref<siconos::algebra::SiconosMatrix> vertices)
    : SiconosShape(), _vertices(std::make_shared<siconos::algebra::SiconosMatrix>(vertices)) {
  _indexes = std::make_shared<std::vector<unsigned int>>(indexes);
  if (!_indexes || (_indexes->size() % 3) != 0)
    THROW_EXCEPTION("Mesh indexes size must be divisible by 3.");
  if (!_vertices || _vertices->rows() != 3)
    THROW_EXCEPTION("Mesh vertices matrix must have 3 columns.");
}

siconos::collision::SiconosBox2d::SiconosBox2d(double width, double height)
    : SiconosShape(), _dimensions(std::make_shared<siconos::algebra::SiconosVector>(2)) {
  (*_dimensions)(0) = width;
  (*_dimensions)(1) = height;
}

void siconos::collision::SiconosBox2d::setDimensions(double width, double height) {
  (*_dimensions)(0) = width;
  (*_dimensions)(1) = height;
  _version++;
}

void siconos::collision::SiconosBox2d::setDimensions(
    const siconos::algebra::SiconosVector& dim) {
  (*_dimensions)(0) = dim(0);
  (*_dimensions)(1) = dim(1);
  _version++;
}

siconos::collision::SiconosConvexHull2d::SiconosConvexHull2d(
    std::shared_ptr<siconos::algebra::SiconosMatrix> vertices)
    : SiconosShape(), _vertices(vertices) {
  if (_vertices && _vertices->cols() != 2)
    THROW_EXCEPTION("Convex hull vertices matrix must have 2 columns in 2d.");
}

siconos::collision::SiconosConvexHull2d::SiconosConvexHull2d(
    Eigen::Ref<siconos::algebra::SiconosMatrix> vertices)
    : SiconosShape(), _vertices(std::make_shared<siconos::algebra::SiconosMatrix>(vertices)) {
  if (_vertices && _vertices->cols() != 2)
    THROW_EXCEPTION("Convex hull vertices matrix must have 2 columns in 2d.");
}

void siconos::collision::SiconosPlane::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosSphere::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosBox::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosCylinder::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosCone::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosCapsule::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosConvexHull::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosMesh::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosHeightMap::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosDisk::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosBox2d::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
void siconos::collision::SiconosConvexHull2d::accept(
    std::shared_ptr<siconos::collision::internal::ShapeVisitor> tourist) {
  tourist->visit(this->shared_from_this());
}
