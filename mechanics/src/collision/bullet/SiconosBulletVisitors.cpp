// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SiconosBulletVisitors.hpp"

#include "RigidBody2dDS.hpp"
#include "RigidBodyDS.hpp"
#include "SiconosBulletCollisionManager_impl.hpp"
#include "SiconosShape.hpp"  // for SiconosPlane, SiconosBox and so on

void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosPlane> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosSphere> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosBox> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosCylinder> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosCone> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosCapsule> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosConvexHull> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosMesh> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<SiconosHeightMap> shape) {
  auto rbds = std::static_pointer_cast<RigidBodyDS>(ds);
  impl->createCollisionObject(base, rbds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosDisk> shape) {
  auto rb2dds = std::static_pointer_cast<RigidBody2dDS>(ds);
  impl->createCollisionObject(base, rb2dds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosBox2d> shape) {
  auto rb2dds = std::static_pointer_cast<RigidBody2dDS>(ds);
  impl->createCollisionObject(base, rb2dds, shape, contactor, staticBody);
}
void siconos::collision::bullet::internal::CreateCollisionObjectShapeVisitor::visit(
    std::shared_ptr<siconos::collision::SiconosConvexHull2d> shape) {
  auto rb2dds = std::static_pointer_cast<RigidBody2dDS>(ds);
  impl->createCollisionObject(base, rb2dds, shape, contactor, staticBody);
}
