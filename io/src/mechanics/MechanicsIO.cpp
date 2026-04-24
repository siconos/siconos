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

#include "MechanicsIO.hpp"

#include "SiconosConfig.h"
#include "Tools.hpp"
#ifdef SICONOS_HAS_BULLET
#include "Bullet2d3DR.hpp"
#include "Bullet2dR.hpp"
#include "Bullet5DR.hpp"
#include "BulletR.hpp"
#else
#include "NewtonEuler3DR.hpp"
#include "NewtonEuler5DR.hpp"
#include "SpaceFilter.hpp"
#endif
#include <concepts>

#include "BlockVector.hpp"
#include "BodyShapeRecord.hpp"
#include "Circle.hpp"
#include "CircleCircleR.hpp"
#include "Contact2d3DR.hpp"
#include "Contact2dR.hpp"
#include "Contact5DR.hpp"
#include "ContactR.hpp"
#include "Disk.hpp"
#include "DiskDiskR.hpp"
#include "DiskPlanR.hpp"
#include "DynamicalSystem.hpp"
#include "FremondImpactFrictionNSL.hpp"
#include "Interaction.hpp"
#include "Lagrangian2d2DR.hpp"
#include "Lagrangian2d3DR.hpp"
#include "LagrangianDS.hpp"
#include "LagrangianSparseDS.hpp"
#include "NSLVisitor.hpp"
#include "NewtonEuler1DR.hpp"
#include "NewtonEuler3DR.hpp"
#include "NewtonEuler5DR.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactRollingFrictionNSL.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "PivotJointR.hpp"
#include "PrismaticJointR.hpp"
#include "Question.hpp"
#include "RigidBody2dDS.hpp"
#include "RigidBodyDS.hpp"
#include "SiconosJoints.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "SimulationGraphs.hpp"
#include "StaticBody.hpp"
#include "Topology.hpp"
// /* all the visitable classes must have been included at this point */
#include "VisitorMaker.hpp"
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

// -------------------
//  Visitor stuff to collect all states, velocities, ... of the nsds and put them in a matrix
// ------------------

struct siconos::io::GetPosition : public siconos::modeling::dynamical_systems::Visitor {
  std::shared_ptr<siconos::algebra::MapVectorType> result;
  void setMap(Eigen::Ref<siconos::algebra::SiconosVector> buffer) override {
    result = std::make_shared<siconos::algebra::MapVectorType>(buffer.data(), buffer.size());
  }

  template <typename T>
  void operator()(const T& ds) {
    (*result)(0) = static_cast<siconos::algebra::SiconosVector::Scalar>(ds.number());
    result->segment(1, ds.q_read().size()) = ds.q_read();
  }
};

struct siconos::io::GetVelocity : public siconos::modeling::dynamical_systems::Visitor {
  std::shared_ptr<siconos::algebra::MapVectorType> result;
  void setMap(Eigen::Ref<siconos::algebra::SiconosVector> buffer) override {
    result = std::make_shared<siconos::algebra::MapVectorType>(buffer.data(), buffer.size());
  }
  template <typename T>
  void operator()(const T& ds) {
    (*result)(0) = static_cast<siconos::algebra::SiconosVector::Scalar>(ds.number());
    result->segment(1, ds.velocity_read().size()) = ds.velocity_read();
  }
};

struct siconos::io::ForMu : public siconos::modeling::nonsmooth_laws::Question<double> {
  using siconos::modeling::nonsmooth_laws::Visitor::visit;

  void visit(const siconos::modeling::NewtonImpactFrictionNSL& nsl) override {
    answer = nsl.mu();
  }
  void visit(const siconos::modeling::FremondImpactFrictionNSL& nsl) override {
    answer = nsl.mu();
  }
  void visit(const siconos::modeling::NewtonImpactRollingFrictionNSL& nsl) override {
    answer = nsl.mu();
  }
  void visit(const siconos::modeling::NewtonImpactNSL& nsl) override { answer = 0.; }
};

struct siconos::io::ForE : public siconos::modeling::nonsmooth_laws::Question<double> {
  using siconos::modeling::nonsmooth_laws::Visitor::visit;

  void visit(const siconos::modeling::NewtonImpactFrictionNSL& nsl) override {
    answer = nsl.en();
  }
  void visit(const siconos::modeling::FremondImpactFrictionNSL& nsl) override {
    answer = nsl.en();
  }
  void visit(const siconos::modeling::NewtonImpactRollingFrictionNSL& nsl) override {
    answer = nsl.en();
  }
  void visit(const siconos::modeling::NewtonImpactNSL& nsl) override { answer = 0.; }
};

/* Get contact informations */
/* default: a visitor that do nothing */
struct siconos::io::ContactPointVisitor : public siconos::modeling::relations::Visitor {
  std::shared_ptr<siconos::modeling::Interaction> inter{nullptr};
  siconos::algebra::SiconosVector answer;

  template <typename T>
  void operator()(const T& rel) {}
};

/* then specializations : */
template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::modeling::NewtonEuler3DR& rel) {
  const auto& posa = rel.pc1();
  const auto& posb = rel.pc2();
  const auto& nc = rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));
  DEBUG_PRINTF("posa(2)=%g\n", posa(2));

  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  siconos::algebra::SiconosVector cf{rel.H_NE_prod_T().cols()};
  cf.noalias() = rel.H_NE_prod_T().transpose() * *inter->lambda(1);
  answer.resize(23);

  answer(0) = mu;
  answer(1) = posa(0);
  answer(2) = posa(1);
  answer(3) = posa(2);
  answer(4) = posb(0);
  answer(5) = posb(1);
  answer(6) = posb(2);
  answer(7) = nc(0);
  answer(8) = nc(1);
  answer(9) = nc(2);
  answer(10) = cf(0);
  answer(11) = cf(1);
  answer(12) = cf(2);
  answer(13) = (*inter->y(0))(0);
  answer(14) = (*inter->y(0))(1);
  answer(15) = (*inter->y(0))(2);
  answer(16) = (*inter->y(1))(0);
  answer(17) = (*inter->y(1))(1);
  answer(18) = (*inter->y(1))(2);
  answer(19) = (*inter->lambda(1))(0);
  answer(20) = (*inter->lambda(1))(1);
  answer(21) = (*inter->lambda(1))(2);
  answer(22) = id;
}

/* then specializations : */
template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::modeling::NewtonEuler5DR& rel) {
  const auto& posa = rel.pc1();
  const auto& posb = rel.pc2();
  const auto& nc = rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));
  DEBUG_PRINTF("posa(2)=%g\n", posa(2));

  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  siconos::algebra::SiconosVector cf{rel.H_NE_prod_T().cols()};
  cf.noalias() = rel.H_NE_prod_T().transpose() * *inter->lambda(1);
  answer.resize(23);

  answer(0) = mu;
  answer(1) = posa(0);
  answer(2) = posa(1);
  answer(3) = posa(2);
  answer(4) = posb(0);
  answer(5) = posb(1);
  answer(6) = posb(2);
  answer(7) = nc(0);
  answer(8) = nc(1);
  answer(9) = nc(2);
  answer(10) = cf(0);
  answer(11) = cf(1);
  answer(12) = cf(2);
  answer(13) = (*inter->y(0))(0);
  answer(14) = (*inter->y(0))(1);
  answer(15) = (*inter->y(0))(2);
  answer(16) = (*inter->y(1))(0);
  answer(17) = (*inter->y(1))(1);
  answer(18) = (*inter->y(1))(2);
  answer(19) = (*inter->lambda(1))(0);
  answer(20) = (*inter->lambda(1))(1);
  answer(21) = (*inter->lambda(1))(2);
  answer(22) = id;
}

template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::collision::native::bodies::DiskDiskR& rel) {
  const auto& ds_vars = inter->read_dynamical_systems_variables();
  auto& q = *ds_vars[tools::enum_to_index(modeling::LagrangianR::ds_var::q0)];

  auto x1 = q(0);
  auto y1 = q(1);
  const auto r1 = rel.getRadius1();
  auto x2 = q(3);
  auto y2 = q(4);
  const auto r2 = rel.getRadius2();
  auto d = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

  auto ncx = d > 0 ? (x2 - x1) / d : 0;
  auto ncy = d > 0 ? (y2 - y1) / d : 0;

  auto cpax = x1 + ncx * r1;
  auto cpay = y1 + ncy * r1;

  auto cpbx = x2 - ncx * r2;
  auto cpby = y2 - ncy * r2;

  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  siconos::algebra::SiconosVector cf{rel.jacobianhOver_q().cols()};
  cf.noalias() = rel.jacobianhOver_q().transpose() * *inter->lambda(1);
  answer.resize(16);

  answer(0) = mu;
  answer(1) = cpax;
  answer(2) = cpay;

  answer(3) = cpbx;
  answer(4) = cpby;

  answer(5) = ncx;
  answer(6) = ncy;

  answer(7) = cf(0);
  answer(8) = cf(1);

  answer(9) = (*inter->y(0))(0);
  answer(10) = (*inter->y(0))(1);

  answer(11) = (*inter->y(1))(0);
  answer(12) = (*inter->y(1))(1);

  answer(13) = (*inter->lambda(1))(0);
  answer(14) = (*inter->lambda(1))(1);

  answer(15) = id;
}

// CircleCircleR should be named DiskCircleR
template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::collision::native::bodies::CircleCircleR& rel) {
  const auto& ds_vars = inter->read_dynamical_systems_variables();
  auto& q = *ds_vars[tools::enum_to_index(modeling::LagrangianR::ds_var::q0)];

  auto x1 = q(0);
  auto y1 = q(1);
  const auto r1 = rel.getRadius1();
  auto x2 = q(3);
  auto y2 = q(4);
  const auto r2 = rel.getRadius2();
  auto d = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

  auto ncx = d > 0 ? (x2 - x1) / d : 0;
  auto ncy = d > 0 ? (y2 - y1) / d : 0;

  double cpax, cpay, cpbx, cpby;
  if (r1 < r2)  // disk1 inside circle2
  {
    cpax = x1 - ncx * r1;
    cpay = y1 - ncy * r1;

    cpbx = x2 - ncx * r2;
    cpby = y2 - ncy * r2;
  } else  // disk2 inside circle1
  {
    cpbx = x2 + ncx * r2;
    cpby = y2 + ncx * r2;

    cpax = x1 + ncx * r1;
    cpay = y1 + ncy * r1;
  }
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  siconos::algebra::SiconosVector cf{rel.jacobianhOver_q().cols()};
  cf.noalias() = rel.jacobianhOver_q().transpose() * *inter->lambda(1);
  answer.resize(16);

  answer(0) = mu;
  answer(1) = cpax;
  answer(2) = cpay;

  answer(3) = cpbx;
  answer(4) = cpby;

  answer(5) = ncx;
  answer(6) = ncy;

  answer(7) = cf(0);
  answer(8) = cf(1);

  answer(9) = (*inter->y(0))(0);
  answer(10) = (*inter->y(0))(1);

  answer(11) = (*inter->y(1))(0);
  answer(12) = (*inter->y(1))(1);

  answer(13) = (*inter->lambda(1))(0);
  answer(14) = (*inter->lambda(1))(1);

  answer(15) = id;
}

template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::collision::native::bodies::DiskPlanR& rel) {
  const auto& ds_vars = inter->read_dynamical_systems_variables();
  const auto& q0 = *ds_vars[tools::enum_to_index(modeling::LagrangianR::ds_var::q0)];

  auto x1 = q0(0);
  auto y1 = q0(1);
  auto r1 = rel.getRadius();
  auto A = rel.getA();
  auto B = rel.getB();
  auto C = rel.getC();

  auto x2 = -(A * C - B * B * x1 + A * B * y1) / (A * A + B * B);
  auto y2 = -(B * C - A * A * y1 + A * B * x1) / (A * A + B * B);

  auto d = (*(inter->y()[1]))(0) + r1;

  auto ncx = d > 0 ? (x2 - x1) / d : 0;
  auto ncy = d > 0 ? (y2 - y1) / d : 0;

  auto cpax = x1 + ncx * r1;
  auto cpay = y1 + ncy * r1;

  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  siconos::algebra::SiconosVector cf{rel.jacobianhOver_q().cols()};
  cf.noalias() = rel.jacobianhOver_q().transpose() * *inter->lambda(1);
  answer.resize(16);

  answer(0) = mu;
  answer(1) = cpax;
  answer(2) = cpay;

  answer(3) = x2;
  answer(4) = y2;

  answer(5) = ncx;
  answer(6) = ncy;

  answer(7) = cf(0);
  answer(8) = cf(1);

  answer(9) = (*inter->y(0))(0);
  answer(10) = (*inter->y(0))(1);

  answer(11) = (*inter->y(1))(0);
  answer(12) = (*inter->y(1))(1);

  answer(13) = (*inter->lambda(1))(0);
  answer(14) = (*inter->lambda(1))(1);

  answer(15) = id;
}

template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::modeling::Lagrangian2d2DR& rel) {
  const auto& posa = rel.pc1();
  const auto& posb = rel.pc2();
  const auto& nc = rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));

  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  siconos::algebra::SiconosVector cf{rel.jacobianhOver_q().cols()};
  cf.noalias() = rel.jacobianhOver_q().transpose() * *inter->lambda(1);

  answer.resize(16);

  answer(0) = mu;
  answer(1) = posa(0);
  answer(2) = posa(1);

  answer(3) = posb(0);
  answer(4) = posb(1);

  answer(5) = nc(0);
  answer(6) = nc(1);

  answer(7) = cf(0);
  answer(8) = cf(1);

  answer(9) = (*inter->y(0))(0);
  answer(10) = (*inter->y(0))(1);

  answer(11) = (*inter->y(1))(0);
  answer(12) = (*inter->y(1))(1);

  answer(13) = (*inter->lambda(1))(0);
  answer(14) = (*inter->lambda(1))(1);

  answer(15) = id;
};

template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::modeling::Lagrangian2d3DR& rel) {
  const auto& posa = rel.pc1();
  const auto& posb = rel.pc2();
  const auto& nc = rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));

  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  siconos::algebra::SiconosVector cf{rel.jacobianhOver_q().cols()};
  cf.noalias() = rel.jacobianhOver_q().transpose() * *inter->lambda(1);

  answer.resize(16);

  answer(0) = mu;
  answer(1) = posa(0);
  answer(2) = posa(1);

  answer(3) = posb(0);
  answer(4) = posb(1);

  answer(5) = nc(0);
  answer(6) = nc(1);

  answer(7) = cf(0);
  answer(8) = cf(1);

  answer(9) = (*inter->y(0))(0);
  answer(10) = (*inter->y(0))(1);

  answer(11) = (*inter->y(1))(0);
  answer(12) = (*inter->y(1))(1);

  answer(13) = (*inter->lambda(1))(0);
  answer(14) = (*inter->lambda(1))(1);

  answer(15) = id;
};

struct siconos::io::ContactPointDomainVisitor : public siconos::modeling::relations::Visitor {
  std::shared_ptr<siconos::modeling::Interaction> inter{nullptr};
  siconos::algebra::SiconosVector answer;

  template <typename T>
  void operator()(const T& rel) {}
};

template <>
void siconos::io::ContactPointDomainVisitor::operator()(
    const siconos::modeling::NewtonEuler3DR& rel) {
  answer.resize(2);

  /*
   * TODO: contact point domain coloring (e.g. based on broadphase).
   * currently, domain = (x>0):1?0
   */
  answer(0) = rel.pc1()(0) > 0;

  answer(1) = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
}

siconos::algebra::SiconosMatrix siconos::io::MechanicsIO::domains(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) const {
  if (nsds.topology()->numberOfIndexSet() < 1)
    return siconos::algebra::SiconosMatrix{};  // empty matrix

  auto& graph = *nsds.topology()->indexSet(1);
  siconos::algebra::SiconosMatrix result =
      siconos::algebra::SiconosMatrix::Zero(graph.vertices_number(), 2);
  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  siconos::algebra::Index current_row;
  for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend;
       ++vi, ++current_row) {
    DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

    using DomainInspector = siconos::internal::RelationVisitor<
        siconos::internal::Classes<siconos::modeling::NewtonEuler1DR,
                                   siconos::modeling::NewtonEuler3DR,
                                   siconos::joints::PrismaticJointR,
                                   siconos::joints::KneeJointR, siconos::joints::PivotJointR>,
        ContactPointDomainVisitor>::Make;

    DomainInspector inspector;
    inspector.inter = graph.bundle(*vi);
    graph.bundle(*vi)->relation()->accept(inspector);
    const auto& data = inspector.answer;
    if (data.size() == 2) result.row(current_row) = data;
  }
  return result;  // RVO
}

template <typename T, typename G>
siconos::algebra::SiconosMatrix siconos::io::MechanicsIO::visitAllVerticesForVector(
    const G& graph) const {
  // Temp. Iterate through the graph to count the number of DS and
  // the max size of these systems
  typename G::VIterator vi, viend;
  siconos::algebra::Index current_col;

  siconos::algebra::Index max_size = 0;
  for (current_col = 0, std::tie(vi, viend) = graph.vertices(); vi != viend;
       ++vi, ++current_col) {
    auto ds = std::dynamic_pointer_cast<siconos::modeling::SecondOrderDS>(graph.bundle(*vi));
    max_size = std::max(max_size, (ds->q_read()).size());
  }

  siconos::algebra::SiconosMatrix results =
      siconos::algebra::SiconosMatrix::Zero(max_size + 1, current_col);

  for (current_col = 0, std::tie(vi, viend) = graph.vertices(); vi != viend;
       ++vi, ++current_col) {
    T getter;
    getter.setMap(results.col(current_col));
    graph.bundle(*vi)->accept(getter);
  }
  return results;  // RVO, no copy
}

template <typename T, typename G>
siconos::algebra::SiconosVector siconos::io::MechanicsIO::visitAllVerticesForDouble(
    const G& graph) const {
  siconos::algebra::SiconosVector result{graph.vertices_number()};

  typename G::VIterator vi, viend;
  siconos::algebra::Index current_row;
  for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend;
       ++vi, ++current_row) {
    T getter;
    graph.bundle(*vi)->accept(getter);
    result(current_row) = getter.result;
  }
  return result;
}

siconos::algebra::SiconosMatrix siconos::io::MechanicsIO::positions(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) const {
  using Getter = siconos::internal::DSVisitor<
      siconos::internal::Classes<siconos::modeling::LagrangianDS,
                                 siconos::modeling::LagrangianSparseDS,
                                 siconos::modeling::NewtonEulerDS>,
      GetPosition>::Make;

  return visitAllVerticesForVector<Getter>(*(nsds.topology()->dSG(0)));
};

siconos::algebra::SiconosMatrix siconos::io::MechanicsIO::velocities(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) const {
  using Getter = siconos::internal::DSVisitor<
      siconos::internal::Classes<siconos::modeling::LagrangianDS,
                                 siconos::modeling::LagrangianSparseDS,
                                 siconos::modeling::NewtonEulerDS>,
      GetVelocity>::Make;

  return visitAllVerticesForVector<Getter>(*nsds.topology()->dSG(0));
}

siconos::algebra::SiconosMatrix siconos::io::MechanicsIO::contactPoints(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds, unsigned int index_set) const {
  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  if (nsds.topology()->numberOfIndexSet() < 1)
    return siconos::algebra::SiconosMatrix{};  // RVO, 0-sized matrix

  // if (nsds.topology()->numberOfIndexSet() > 0) {
  auto& graph = *nsds.topology()->indexSet(index_set);
  siconos::algebra::Index current_row;
  siconos::algebra::SiconosMatrix result =
      siconos::algebra::SiconosMatrix::Zero(graph.vertices_number(), 25);

  int data_size = 0;
  for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend; ++vi) {
    DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

    /* create a visitor for specified classes */
    using ContactPointInspector = siconos::internal::RelationVisitor<
        siconos::internal::Classes<
            siconos::modeling::NewtonEuler1DR, siconos::modeling::NewtonEuler3DR,
            siconos::modeling::NewtonEuler5DR, siconos::modeling::Lagrangian2d2DR,
            siconos::modeling::Lagrangian2d3DR,
            siconos::collision::native::bodies::CircleCircleR,
            siconos::collision::native::bodies::DiskDiskR,
            siconos::collision::native::bodies::DiskPlanR>,
        ContactPointVisitor>::Make;

    ContactPointInspector inspector;
    inspector.inter = graph.bundle(*vi);
    graph.bundle(*vi)->relation()->accept(inspector);
    siconos::algebra::SiconosVector& data = inspector.answer;
    data_size = data.size();

    if (data_size == 0) {
      // Nothing is done since the relation does not appear as a relation
      // related to a contact points (perhaps a joint)
    } else {
      // We add at the end the number of ds1 and ds2
      data.conservativeResize(data_size + 2);
      DEBUG_EXPR(siconos::algebra::print(data););
      auto& ds1 = *graph.properties(*vi).source;
      auto& ds2 = *graph.properties(*vi).target;
      data(data_size) = static_cast<siconos::algebra::SiconosVector::Scalar>(ds1.number());
      data(data_size + 1) = static_cast<siconos::algebra::SiconosVector::Scalar>(ds2.number());
      DEBUG_EXPR(siconos::algebra::print(data););
      if (result.cols() != data.size()) {
        result.conservativeResize(siconos::algebra::to_index(graph.vertices_number()),
                                  data.size());
      }
      result.row(current_row++) = data;
      data_size += 2;
    }
  }
  result.conservativeResize(current_row, data_size);
  DEBUG_EXPR(siconos::algebra::print(result));
  return result;  // RVO
}

/* Get contact informations */
/* default: a visitor that do nothing */
struct ContactInfoVisitor : public siconos::modeling::relations::Visitor {
  std::shared_ptr<siconos::modeling::Interaction> inter{nullptr};
  // std::vector<int> answer; better with a vector of int
  siconos::algebra::SiconosVector answer;

  template <typename T>
  void operator()(const T& rel) {}
};

/* then specializations : */
template <>
void ContactInfoVisitor::operator()(const siconos::modeling::NewtonEuler3DR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(10);
  // answer[0]= id;
  // answer[1]= 0; // reserve for ds1.number
  // answer[2]= 0; // reserve for ds2.number
  answer(0) = id;
  answer(1) = 0;
  answer(2) = 0;
}

template <>
void ContactInfoVisitor::operator()(const siconos::collision::ContactR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(4);
  answer(0) = id;
  answer(1) = 0;
  answer(2) = 0;
  if (rel.bodyShapeRecordB->staticBody) {
    answer(3) = rel.bodyShapeRecordB->staticBody->number;
  } else
    answer(3) = 0;
}

template <>
void ContactInfoVisitor::operator()(const siconos::collision::Contact5DR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(4);
  answer(0) = id;
  answer(1) = 0;
  answer(2) = 0;
  if (rel.bodyShapeRecordB->staticBody) {
    answer(3) = rel.bodyShapeRecordB->staticBody->number;
  } else
    answer(3) = 0;
}

template <>
void ContactInfoVisitor::operator()(const siconos::collision::Contact2dR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(4);
  answer(0) = id;
  answer(1) = 0;
  answer(2) = 0;
  if (rel.bodyShapeRecordB->staticBody) {
    answer(3) = rel.bodyShapeRecordB->staticBody->number;
  } else
    answer(3) = 0;
}

template <>
void ContactInfoVisitor::operator()(const siconos::collision::Contact2d3DR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(4);
  answer(0) = id;
  answer(1) = 0;
  answer(2) = 0;
  if (rel.bodyShapeRecordB->staticBody) {
    answer(3) = rel.bodyShapeRecordB->staticBody->number;
  } else
    answer(3) = 0;
}

std::optional<siconos::algebra::SiconosMatrix> siconos::io::MechanicsIO::contactInfo(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds, unsigned int index_set) const {
  DEBUG_BEGIN("MechanicsIO::contactInfo");

  siconos::graphs::InteractionsGraph::VIterator vi, viend;

  if (!(nsds.topology()->numberOfIndexSet() > 0)) return std::nullopt;

  auto& graph = *nsds.topology()->indexSet(index_set);
  siconos::algebra::Index current_row;
  siconos::algebra::SiconosMatrix result =
      siconos::algebra::SiconosMatrix::Zero(graph.vertices_number(), 4);

  int data_size = 0;
  for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend; ++vi) {
    DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

    /* create a visitor for specified classes */
    using ContactInfoInspector = siconos::internal::RelationVisitor<
        siconos::internal::Classes<
            siconos::modeling::NewtonEuler3DR, siconos::collision::ContactR,
            siconos::collision::Contact5DR, siconos::collision::Contact2dR,
            siconos::collision::Contact2d3DR>,
        ContactInfoVisitor>::Make;

    ContactInfoInspector inspector;
    inspector.inter = graph.bundle(*vi);
    graph.bundle(*vi)->relation()->accept(inspector);
    auto& data = inspector.answer;
    data_size = data.size();

    if (data_size == 0) {
      // Nothing is done since the relation does not appear as a relation
      // related to a contact points (perhaps a joint)
    } else {
      // We add at the end the number of ds1 and ds2
      DEBUG_EXPR(siconos::algebra::print(data););
      auto& ds1 = *graph.properties(*vi).source;
      auto& ds2 = *graph.properties(*vi).target;
      data(1) = static_cast<siconos::algebra::SiconosVector::Scalar>(ds1.number());
      data(2) = static_cast<siconos::algebra::SiconosVector::Scalar>(ds2.number());
    }
    if (result.cols() != data.size()) {
      result.conservativeResize(siconos::algebra::to_index(graph.vertices_number()),
                                data.size());
    }
    result.row(current_row++) = data;
  }
  // result.resize(current_row, data_size);
  // Why do we need a resize? If result is too small, the line 748 will fail anyway.
  // Note FP: warning resize with eigen may lead to reset
  DEBUG_EXPR(siconos::algebra::print(result););
  return result;  // RVO, no copy

  DEBUG_END(" MechanicsIO::contactInfo");
}

namespace siconos::io {
/* Get contact work information */
/* default: a visitor that do nothing */

struct ContactContactWorkVisitor : public siconos::modeling::relations::Visitor {
  std::shared_ptr<siconos::modeling::Interaction> inter{nullptr};
  // std::vector<int> answer; better with a vector of int
  siconos::algebra::SiconosVector answer;
  double omega;
  double tol;
  template <typename T>
  void operator()(const T& rel) {}
};

/* then specializations : */
template <>
void ContactContactWorkVisitor::operator()(const siconos::modeling::NewtonEuler3DR& rel) {
  answer.resize(7);
}

static void compute_contact_work_and_status(
    std::shared_ptr<siconos::modeling::Interaction> inter, double theta, double tol,
    siconos::algebra::SiconosVector& answer) {
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  auto e = siconos::modeling::nonsmooth_laws::ask<ForE>(*inter->nonSmoothLaw());

  // Compute normal contact work
  auto vn_minus = (inter->y_k(1))(0);
  auto vn_plus = (*inter->y(1))(0);
  auto pn = (*inter->lambda(1))(0);

  auto normal_contact_work = 0.5 * (vn_minus + vn_plus) * pn;
  answer(1) = normal_contact_work;

  double vn_average_theta = theta * vn_plus + (1. - theta) * vn_minus;
  double normal_contact_work_theta = vn_average_theta * pn;
  answer(3) = normal_contact_work_theta;

  // Compute tangent contact work of impulse
  auto vt_1_minus = (inter->y_k(1))(1);
  auto vt_2_minus = (inter->y_k(1))(2);
  auto vt_1_plus = (*inter->y(1))(1);
  auto vt_2_plus = (*inter->y(1))(2);

  double vt_1_average = 1 / 2. * (vt_1_plus + vt_1_minus);
  double vt_2_average = 1 / 2. * (vt_2_plus + vt_2_minus);

  double pt_1 = (*inter->lambda(1))(1);
  double pt_2 = (*inter->lambda(1))(2);

  double tangent_contact_work = vt_1_average * pt_1 + vt_2_average * pt_2;
  answer(2) = tangent_contact_work;

  double vt_1_average_theta = theta * vt_1_plus + (1. - theta) * vt_1_minus;
  double vt_2_average_theta = theta * vt_2_plus + (1. - theta) * vt_2_minus;

  double tangent_contact_work_theta = vt_1_average_theta * pt_1 + vt_2_average_theta * pt_2;
  answer(4) = tangent_contact_work_theta;

  // // Compute directly work dissipated by friction impulse
  // double norm_vt_average = sqrt(vt_1_average * vt_1_average + vt_2_average *
  // vt_2_average); double friction_dissipation = mu * norm_vt_average * pn; answer(3) =
  // friction_dissipation;

  /* Compute contact status
   * Warning the status are consistent for the sticking and sliding
   * only with fully implicit discretization o NewtonImpact law
   * and not wih Fremond impact law
   */

  double norm_pt = sqrt(pt_1 * pt_1 + pt_2 * pt_2);
  double norm_vt_plus = sqrt(vt_1_plus * vt_1_plus + vt_2_plus * vt_2_plus);
  // double norm_vt_minus = sqrt(vt_1_minus * vt_1_minus + vt_2_minus * vt_2_minus);
  if ((pn < tol) and (vn_plus + e * vn_minus > tol))
    answer(5) = 0;  // take-off = 0
  else if ((pn > tol) and (vn_plus + e * vn_minus < tol)) {
    if ((norm_pt - mu * pn > tol)) {
      // std::cout << "WARNING: the impulse is outside the Coulomb cone  " << std::endl;
      answer(5) = -3;  // outside the cone = -3
    } else if ((norm_pt - mu * pn < -tol)) {
      // std::cout << "the impulse is in the *interior* of  the Coulomb cone  " << std::endl;
      // std::cout << "norm_vt_plus  " << norm_vt_plus << std::endl;
      if (norm_vt_plus > tol) {
        // std::cout << "WARNING: but the norm of vt is not zero  " << std::endl;
        answer(5) = -2;  // sticking with a non zero slifing velocity = -2
      }
      answer(5) = 1;  // sticking = 1
    } else {
      // std::cout << "the impulse is on the *boundary* of the Coulomb cone  " << std::endl;
      // std::cout << "norm_vt_plus  " << norm_vt_plus << std::endl;
      answer(5) = 2;  // sliding = 2
    }
  } else
    answer(5) = -1;  // undetermined = -1

  if ((pn > tol) and (vn_minus > tol)) {
    // std::cout << "WARNING: we apply the impact law of positive velocity " << std::endl;
    // std::cout << "pn " << pn << " vn minus " << vn_minus << " vn plus " << vn_plus
    // 		<< " normal_contact_work " << normal_contact_work
    // 		<< " -e * vn_minus   " << -e*vn_minus
    // 		<< std::endl;
    answer(6) = normal_contact_work;
  }
  //   double id = inter->number();
  //   std::cout << "\nid "<< id << " 3D" << std::endl;
  //   std::cout << " e "<< e  << " mu "<< mu << std::endl;
  //   std::cout << " tol "<< tol << " theta " << theta << std::endl;
  //   std::cout << "vn_plus "<< vn_plus << std::endl;
  //   std::cout << "vn_minus "<< vn_minus << std::endl;
  //   std::cout << "pn "<< pn << std::endl;
  //   std::cout << "normal_contact_work  "<< normal_contact_work  << std::endl;

  //   std::cout << "vt_plus "<< vt_1_plus << " " << vt_2_plus <<  std::endl;
  //   std::cout << "vt_minus "<< vt_1_minus << " " << vt_2_minus <<  std::endl;
  //   std::cout << "pt "<< pt_1 << " "  << pt_2 << std::endl;
  //   std::cout << "tangent_contact_work  "<< tangent_contact_work  << std::endl;

  //   std::cout << "friction_dissipation  "<< friction_dissipation << std::endl;

  //   std::cout << "norm_pt  "<< norm_pt  << std::endl;
  //   std::cout << "norm_pt - mu* pn  "<< norm_pt -mu *pn   << std::endl;
  //   std::cout << "vn_plus + e * vn_minus  " << vn_plus + e * vn_minus   << std::endl;
  //   std::cout << "status   "<<   answer(4) << std::endl;
}

static void compute_contact_work_and_status_2d(
    std::shared_ptr<siconos::modeling::Interaction> inter, double theta, double tol,
    siconos::algebra::SiconosVector& answer) {
  auto mu = siconos::modeling::nonsmooth_laws::ask<ForMu>(*inter->nonSmoothLaw());
  auto e = siconos::modeling::nonsmooth_laws::ask<ForE>(*inter->nonSmoothLaw());

  // Compute normal contact work
  double vn_minus = inter->y_k(1)(0);
  double vn_plus = (*inter->y(1))(0);
  double pn = (*inter->lambda(1))(0);

  double vn_average = 1 / 2. * (vn_plus + vn_minus);
  double normal_contact_work = vn_average * pn;
  answer(1) = normal_contact_work;

  double vn_average_theta = theta * vn_plus + (1. - theta) * vn_minus;
  double normal_contact_work_theta = vn_average_theta * pn;
  answer(3) = normal_contact_work_theta;

  // Compute tangent contact work of impulse
  double vt_1_minus = inter->y_k(1)(1);
  double vt_1_plus = (*inter->y(1))(1);
  double pt_1 = (*inter->lambda(1))(1);

  double vt_1_average = 1 / 2. * (vt_1_plus + vt_1_minus);
  double tangent_contact_work = vt_1_average * pt_1;
  answer(2) = tangent_contact_work;

  double vt_1_average_theta = theta * vt_1_plus + (1. - theta) * vt_1_minus;
  double tangent_contact_work_theta = vt_1_average_theta * pt_1;
  answer(4) = tangent_contact_work_theta;

  // // Compute directly work dissipated by friction impulse
  // double norm_vt_average = sqrt(vt_1_average * vt_1_average);
  // double friction_dissipation = mu * norm_vt_average * pn;
  // answer(3) = friction_dissipation;

  /* Compute contact status
   * Warning the status are consistent for the sticking and sliding
   * only with fully implicit discretization o NewtonImpact law
   * and not wih Fremond impact law
   */

  double norm_pt = sqrt(pt_1 * pt_1);
  double norm_vt_plus = sqrt(vt_1_plus * vt_1_plus);
  // double norm_vt_minus = sqrt(vt_1_minus*vt_1_minus);
  if ((pn < tol) and (vn_plus + e * vn_minus > tol))
    answer(5) = 0;  // take-off = 0
  else if ((pn > tol) and (vn_plus + e * vn_minus < tol)) {
    if ((norm_pt - mu * pn > tol)) {
      // std::cout << "WARNING: the impulse is outside the Coulomb cone  " << std::endl;
      answer(5) = -3;  // outside the cone = -3
    } else if ((norm_pt - mu * pn < -tol)) {
      // std::cout << "the impulse is in the *interior* of  the Coulomb cone  " << std::endl;
      // std::cout << "norm_vt_plus  " << norm_vt_plus << std::endl;
      if (norm_vt_plus > tol) {
        // std::cout << "WARNING: but the norm of vt is not zero  " << std::endl;
        answer(5) = -2;  // sticking with a non zero slifing velocity = -2
      }
      // ?? answer(5) = -2 is always overwritten
      answer(5) = 1;  // sticking = 1
    } else {
      // std::cout << "the impulse is on the *boundary* of the Coulomb cone  " << std::endl;
      // std::cout << "norm_vt_plus  " << norm_vt_plus << std::endl;
      answer(5) = 2;  // sliding = 2
    }
  } else
    answer(5) = -1;  // undetermined = -1
  if ((pn > tol) and (vn_minus > tol)) {
    // std::cout << "WARNING: we apply the impact law of positive velocity " << std::endl;
    // std::cout << "pn " << pn << " vn minus " << vn_minus << " vn plus " << vn_plus
    // 		<< " normal_contact_work " << normal_contact_work
    // 		<< " -e * vn_minus   " << -e*vn_minus
    // 		<< std::endl;
    answer(6) = normal_contact_work;
  }
  // double id = inter->number();
  // std::cout << "\nid "<< id << " 2D" << std::endl;
  // std::cout << " e "<< e  << " mu "<< mu << std::endl;
  // std::cout << " tol "<< tol << " theta " << theta << std::endl;
  // std::cout << "vn_plus "<< vn_plus << std::endl;
  // std::cout << "vn_minus "<< vn_minus << std::endl;
  // std::cout << "pn "<< pn << std::endl;
  // std::cout << "normal_contact_work  "<< normal_contact_work  << std::endl;

  // std::cout << "vt_plus "<< vt_1_plus <<  std::endl;
  // std::cout << "vt_minus "<< vt_1_minus <<  std::endl;
  // std::cout << "pt "<< pt_1 << " "  <<  std::endl;
  // std::cout << "tangent_contact_work  "<< tangent_contact_work  << std::endl;

  // std::cout << "friction_dissipation  "<< friction_dissipation << std::endl;

  // std::cout << "norm_pt  "<< norm_pt  << std::endl;
  // std::cout << "norm_pt - mu* pn  "<< norm_pt -mu *pn   << std::endl;
  // std::cout << "vn_plus + e * vn_minus  " << vn_plus + e * vn_minus   << std::endl;
  // std::cout << "status   "<<   answer(4) << std::endl;
}

}  // namespace siconos::io

template <>
void siconos::io::ContactContactWorkVisitor::operator()(
    const siconos::collision::ContactR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(7);
  answer(0) = id;
  compute_contact_work_and_status(inter, omega, tol, answer);
}

template <>
void siconos::io::ContactContactWorkVisitor::operator()(
    const siconos::collision::Contact5DR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(7);
  answer(0) = id;

  compute_contact_work_and_status(inter, omega, tol, answer);
}

template <>
void siconos::io::ContactContactWorkVisitor::operator()(
    const siconos::collision::Contact2dR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(7);
  answer(0) = id;

  compute_contact_work_and_status_2d(inter, omega, tol, answer);
}

template <>
void siconos::io::ContactContactWorkVisitor::operator()(
    const siconos::collision::Contact2d3DR& rel) {
  auto id = static_cast<siconos::algebra::SiconosVector::Scalar>(inter->number());
  answer.resize(7);
  answer(0) = id;

  compute_contact_work_and_status_2d(inter, omega, tol, answer);
}

siconos::algebra::SiconosMatrix siconos::io::MechanicsIO::contactContactWork(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds, unsigned int index_set,
    double omega, double tol) const {
  DEBUG_BEGIN("SiconosMatrix MechanicsIO::contactContactWork");

  siconos::graphs::InteractionsGraph::VIterator vi, viend;

  if (nsds.topology()->numberOfIndexSet() < 1) return siconos::algebra::SiconosMatrix{};

  auto& graph = *nsds.topology()->indexSet(index_set);

  siconos::algebra::Index current_row;
  siconos::algebra::SiconosMatrix result =
      siconos::algebra::SiconosMatrix::Zero(graph.vertices_number(), 25);

  int data_size = 0;
  for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend; ++vi) {
    DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

    /* create a visitor for specified classes */
    using ContactContactWorkInspector = siconos::internal::RelationVisitor<
        siconos::internal::Classes<
            siconos::modeling::NewtonEuler1DR, siconos::modeling::NewtonEuler3DR,
            siconos::modeling::NewtonEuler5DR, siconos::modeling::Lagrangian2d2DR,
            siconos::modeling::Lagrangian2d3DR, siconos::collision::ContactR,
            siconos::collision::Contact5DR, siconos::collision::Contact2dR,
            siconos::collision::Contact2d3DR,
            siconos::collision::native::bodies::CircleCircleR,
            siconos::collision::native::bodies::DiskDiskR,
            siconos::collision::native::bodies::DiskPlanR>,
        ContactContactWorkVisitor>::Make;
    ContactContactWorkInspector inspector;
    inspector.inter = graph.bundle(*vi);
    inspector.tol = tol;
    inspector.omega = omega;
    graph.bundle(*vi)->relation()->accept(inspector);
    auto& data = inspector.answer;
    data_size = data.size();

    if (data_size == 0) {
      // Nothing is done since the relation does not appear as a relation
      // related to a contact points (perhaps a joint)
    } else {
    }
    if (result.cols() != data.size()) {
      result.conservativeResize(siconos::algebra::to_index(graph.vertices_number()),
                                data.size());
    }
    result.row(current_row++) = data;
  }
  result.conservativeResize(current_row, data_size);
  DEBUG_EXPR(siconos::algebra::print(*result));
  ;
  return result;  // RVO

  DEBUG_END("MechanicsIO::contactContactWork");
}
