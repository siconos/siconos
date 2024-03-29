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
#include "KneeJointR.hpp"
#include "Lagrangian2d2DR.hpp"
#include "Lagrangian2d3DR.hpp"
#include "LagrangianDS.hpp"
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
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "SiconosVisitor.hpp"
#include "SimpleMatrix.hpp"
#include "SimulationGraphs.hpp"
#include "StaticBody.hpp"
#include "Topology.hpp"
// /* all the visitable classes must have been included at this point */
#include "VisitorMaker.hpp"
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

struct siconos::io::GetPosition : public siconos::internal::SiconosVisitor {
  std::shared_ptr<siconos::algebra::SiconosVector> result{nullptr};

  template <typename T>
  void operator()(const T& ds) {
    result = std::make_shared<siconos::algebra::SiconosVector>(1 + ds.q()->size());
    result->setValue(0, ds.number());
    result->setBlock(1, *ds.q());
  }
};

struct siconos::io::GetVelocity : public siconos::internal::SiconosVisitor {
  std::shared_ptr<siconos::algebra::SiconosVector> result{nullptr};

  template <typename T>
  void operator()(const T& ds) {
    result = std::make_shared<siconos::algebra::SiconosVector>(1 + ds.velocity()->size());
    result->setValue(0, ds.number());
    result->setBlock(1, *ds.velocity());
  }
};

struct siconos::io::ForMu : public siconos::internal::Question<double> {
  using SiconosVisitor::visit;

  void visit(const siconos::modeling::NewtonImpactFrictionNSL& nsl) { answer = nsl.mu(); }
  void visit(const siconos::modeling::FremondImpactFrictionNSL& nsl) { answer = nsl.mu(); }
  void visit(const siconos::modeling::NewtonImpactRollingFrictionNSL& nsl) {
    answer = nsl.mu();
  }
  void visit(const siconos::modeling::NewtonImpactNSL& nsl) { answer = 0.; }
};

struct siconos::io::ForE : public siconos::internal::Question<double> {
  using SiconosVisitor::visit;
  void visit(const siconos::modeling::NewtonImpactFrictionNSL& nsl) { answer = nsl.en(); }
  void visit(const siconos::modeling::FremondImpactFrictionNSL& nsl) { answer = nsl.en(); }
  void visit(const siconos::modeling::NewtonImpactRollingFrictionNSL& nsl) {
    answer = nsl.en();
  }
  void visit(const siconos::modeling::NewtonImpactNSL& nsl) { answer = 0.; }
};

/* Get contact informations */
/* default: a visitor that do nothing */
struct siconos::io::ContactPointVisitor : public siconos::internal::SiconosVisitor {
  std::shared_ptr<siconos::modeling::Interaction> inter{nullptr};
  siconos::algebra::SiconosVector answer;

  template <typename T>
  void operator()(const T& rel) {}
};

/* then specializations : */
template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::modeling::NewtonEuler3DR& rel) {
  const auto& posa = *rel.pc1();
  const auto& posb = *rel.pc2();
  const auto& nc = *rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));
  DEBUG_PRINTF("posa(2)=%g\n", posa(2));

  auto id = inter->number();
  auto mu = siconos::internal::ask<ForMu>(*inter->nonSmoothLaw());
  const auto& jachqT = *rel.jachqT();
  siconos::algebra::SiconosVector cf(jachqT.size(1));
  siconos::algebra::prod(*inter->lambda(1), jachqT, cf, true);
  answer.resize(23);

  answer.setValue(0, mu);
  answer.setValue(1, posa(0));
  answer.setValue(2, posa(1));
  answer.setValue(3, posa(2));
  answer.setValue(4, posb(0));
  answer.setValue(5, posb(1));
  answer.setValue(6, posb(2));
  answer.setValue(7, nc(0));
  answer.setValue(8, nc(1));
  answer.setValue(9, nc(2));
  answer.setValue(10, cf(0));
  answer.setValue(11, cf(1));
  answer.setValue(12, cf(2));
  answer.setValue(13, inter->y(0)->getValue(0));
  answer.setValue(14, inter->y(0)->getValue(1));
  answer.setValue(15, inter->y(0)->getValue(2));
  answer.setValue(16, inter->y(1)->getValue(0));
  answer.setValue(17, inter->y(1)->getValue(1));
  answer.setValue(18, inter->y(1)->getValue(2));
  answer.setValue(19, inter->lambda(1)->getValue(0));
  answer.setValue(20, inter->lambda(1)->getValue(1));
  answer.setValue(21, inter->lambda(1)->getValue(2));
  answer.setValue(22, id);
}

/* then specializations : */
template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::modeling::NewtonEuler5DR& rel) {
  const auto& posa = *rel.pc1();
  const auto& posb = *rel.pc2();
  const auto& nc = *rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));
  DEBUG_PRINTF("posa(2)=%g\n", posa(2));

  auto id = inter->number();
  auto mu = siconos::internal::ask<ForMu>(*inter->nonSmoothLaw());
  const auto& jachqT = *rel.jachqT();
  siconos::algebra::SiconosVector cf(jachqT.size(1));
  prod(*inter->lambda(1), jachqT, cf, true);
  answer.resize(23);

  answer.setValue(0, mu);
  answer.setValue(1, posa(0));
  answer.setValue(2, posa(1));
  answer.setValue(3, posa(2));
  answer.setValue(4, posb(0));
  answer.setValue(5, posb(1));
  answer.setValue(6, posb(2));
  answer.setValue(7, nc(0));
  answer.setValue(8, nc(1));
  answer.setValue(9, nc(2));
  answer.setValue(10, cf(0));
  answer.setValue(11, cf(1));
  answer.setValue(12, cf(2));
  answer.setValue(13, inter->y(0)->getValue(0));
  answer.setValue(14, inter->y(0)->getValue(1));
  answer.setValue(15, inter->y(0)->getValue(2));
  answer.setValue(16, inter->y(1)->getValue(0));
  answer.setValue(17, inter->y(1)->getValue(1));
  answer.setValue(18, inter->y(1)->getValue(2));
  answer.setValue(19, inter->lambda(1)->getValue(0));
  answer.setValue(20, inter->lambda(1)->getValue(1));
  answer.setValue(21, inter->lambda(1)->getValue(2));
  answer.setValue(22, id);
}

template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::collision::native::bodies::DiskDiskR& rel) {
  auto& DSlink = inter->linkToDSVariables();
  auto& q = *DSlink[siconos::modeling::LagrangianR::q0];

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

  auto id = inter->number();
  auto mu = siconos::internal::ask<ForMu>(*inter->nonSmoothLaw());
  const auto& jachq = *rel.jachq();
  siconos::algebra::SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);

  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, cpax);
  answer.setValue(2, cpay);

  answer.setValue(3, cpbx);
  answer.setValue(4, cpby);

  answer.setValue(5, ncx);
  answer.setValue(6, ncy);

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9, inter->y(0)->getValue(0));
  answer.setValue(10, inter->y(0)->getValue(1));

  answer.setValue(11, inter->y(1)->getValue(0));
  answer.setValue(12, inter->y(1)->getValue(1));

  answer.setValue(13, inter->lambda(1)->getValue(0));
  answer.setValue(14, inter->lambda(1)->getValue(1));

  answer.setValue(15, id);
}

// CircleCircleR should be named DiskCircleR
template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::collision::native::bodies::CircleCircleR& rel) {
  auto& DSlink = inter->linkToDSVariables();
  auto& q = *DSlink[siconos::modeling::LagrangianR::q0];

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
  auto id = inter->number();
  auto mu = siconos::internal::ask<ForMu>(*inter->nonSmoothLaw());
  const auto& jachq = *rel.jachq();
  siconos::algebra::SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);

  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, cpax);
  answer.setValue(2, cpay);

  answer.setValue(3, cpbx);
  answer.setValue(4, cpby);

  answer.setValue(5, ncx);
  answer.setValue(6, ncy);

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9, inter->y(0)->getValue(0));
  answer.setValue(10, inter->y(0)->getValue(1));

  answer.setValue(11, inter->y(1)->getValue(0));
  answer.setValue(12, inter->y(1)->getValue(1));

  answer.setValue(13, inter->lambda(1)->getValue(0));
  answer.setValue(14, inter->lambda(1)->getValue(1));

  answer.setValue(15, id);
}

template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::collision::native::bodies::DiskPlanR& rel) {
  auto& DSlink = inter->linkToDSVariables();
  const auto& q0 = *DSlink[siconos::modeling::LagrangianR::q0];

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

  auto id = inter->number();
  auto mu = siconos::internal::ask<ForMu>(*inter->nonSmoothLaw());
  const auto& jachq = *rel.jachq();
  siconos::algebra::SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);

  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, cpax);
  answer.setValue(2, cpay);

  answer.setValue(3, x2);
  answer.setValue(4, y2);

  answer.setValue(5, ncx);
  answer.setValue(6, ncy);

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9, inter->y(0)->getValue(0));
  answer.setValue(10, inter->y(0)->getValue(1));

  answer.setValue(11, inter->y(1)->getValue(0));
  answer.setValue(12, inter->y(1)->getValue(1));

  answer.setValue(13, inter->lambda(1)->getValue(0));
  answer.setValue(14, inter->lambda(1)->getValue(1));

  answer.setValue(15, id);
}

template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::modeling::Lagrangian2d2DR& rel) {
  const auto& posa = *rel.pc1();
  const auto& posb = *rel.pc2();
  const auto& nc = *rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));

  auto id = inter->number();
  auto mu = siconos::internal::ask<ForMu>(*inter->nonSmoothLaw());
  const auto& jachq = *rel.jachq();
  siconos::algebra::SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);

  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, posa(0));
  answer.setValue(2, posa(1));

  answer.setValue(3, posb(0));
  answer.setValue(4, posb(1));

  answer.setValue(5, nc(0));
  answer.setValue(6, nc(1));

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9, inter->y(0)->getValue(0));
  answer.setValue(10, inter->y(0)->getValue(1));

  answer.setValue(11, inter->y(1)->getValue(0));
  answer.setValue(12, inter->y(1)->getValue(1));

  answer.setValue(13, inter->lambda(1)->getValue(0));
  answer.setValue(14, inter->lambda(1)->getValue(1));

  answer.setValue(15, id);
};

template <>
void siconos::io::ContactPointVisitor::operator()(
    const siconos::modeling::Lagrangian2d3DR& rel) {
  const auto& posa = *rel.pc1();
  const auto& posb = *rel.pc2();
  const auto& nc = *rel.nc();
  DEBUG_PRINTF("posa(0)=%g\n", posa(0));
  DEBUG_PRINTF("posa(1)=%g\n", posa(1));

  auto id = inter->number();
  auto mu = siconos::internal::ask<ForMu>(*inter->nonSmoothLaw());
  const auto& jachq = *rel.jachq();
  siconos::algebra::SiconosVector cf(jachq.size(1));
  prod(*inter->lambda(1), jachq, cf, true);

  answer.resize(16);

  answer.setValue(0, mu);
  answer.setValue(1, posa(0));
  answer.setValue(2, posa(1));

  answer.setValue(3, posb(0));
  answer.setValue(4, posb(1));

  answer.setValue(5, nc(0));
  answer.setValue(6, nc(1));

  answer.setValue(7, cf(0));
  answer.setValue(8, cf(1));

  answer.setValue(9, inter->y(0)->getValue(0));
  answer.setValue(10, inter->y(0)->getValue(1));

  answer.setValue(11, inter->y(1)->getValue(0));
  answer.setValue(12, inter->y(1)->getValue(1));

  answer.setValue(13, inter->lambda(1)->getValue(0));
  answer.setValue(14, inter->lambda(1)->getValue(1));

  answer.setValue(15, id);
};

struct siconos::io::ContactPointDomainVisitor : public siconos::internal::SiconosVisitor {
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
  answer.setValue(0, rel.pc1()->getValue(0) > 0);

  answer.setValue(1, inter->number());
}

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::io::MechanicsIO::domains(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) const {
  if (nsds.topology()->numberOfIndexSet() > 0) {
    auto& graph = *nsds.topology()->indexSet(1);
    auto result = std::make_shared<siconos::algebra::SimpleMatrix>(graph.vertices_number(), 2);
    siconos::graphs::InteractionsGraph::VIterator vi, viend;
    unsigned int current_row;
    for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend;
         ++vi, ++current_row) {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      using DomainInspector = siconos::internal::Visitor<
          siconos::internal::Classes<
              siconos::modeling::NewtonEuler1DR, siconos::modeling::NewtonEuler3DR,
              siconos::joints::PrismaticJointR, siconos::joints::KneeJointR,
              siconos::joints::PivotJointR>,
          ContactPointDomainVisitor>::Make;

      auto inspector = std::make_shared<DomainInspector>();
      inspector->inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      const auto& data = inspector->answer;
      if (data.size() == 2) result->setRow(current_row, data);
    }
    return result;
  }
  return nullptr;
}

template <typename T, typename G>
std::shared_ptr<siconos::algebra::SimpleMatrix>
siconos::io::MechanicsIO::visitAllVerticesForVector(const G& graph) const {
  auto result = std::make_shared<siconos::algebra::SimpleMatrix>(0, 0);

  typename G::VIterator vi, viend;
  unsigned int current_row;
  for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend;
       ++vi, ++current_row) {
    auto getter = std::make_shared<T>();
    graph.bundle(*vi)->acceptSP(getter);
    const auto& data = *getter->result;
    result->resize(current_row + 1, data.size());
    result->setRow(current_row, data);
  }
  return result;
}

template <typename T, typename G>
std::shared_ptr<siconos::algebra::SiconosVector>
siconos::io::MechanicsIO::visitAllVerticesForDouble(const G& graph) const {
  auto result = std::make_shared<siconos::algebra::SiconosVector>(graph.vertices_number());

  typename G::VIterator vi, viend;
  unsigned int current_row;
  for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend;
       ++vi, ++current_row) {
    auto getter = std::make_shared<T>();
    graph.bundle(*vi)->acceptSP(getter);
    result->setValue(current_row, getter.result);
  }
  return result;
}

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::io::MechanicsIO::positions(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) const {
  using Getter =
      siconos::internal::Visitor<siconos::internal::Classes<siconos::modeling::LagrangianDS,
                                                            siconos::modeling::NewtonEulerDS>,
                                 GetPosition>::Make;

  return visitAllVerticesForVector<Getter>(*(nsds.topology()->dSG(0)));
};

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::io::MechanicsIO::velocities(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) const {
  using Getter =
      siconos::internal::Visitor<siconos::internal::Classes<siconos::modeling::LagrangianDS,
                                                            siconos::modeling::NewtonEulerDS>,
                                 GetVelocity>::Make;

  return visitAllVerticesForVector<Getter>(*nsds.topology()->dSG(0));
}

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::io::MechanicsIO::contactPoints(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds, unsigned int index_set) const {
  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  if (nsds.topology()->numberOfIndexSet() > 0) {
    auto& graph = *nsds.topology()->indexSet(index_set);
    unsigned int current_row;
    auto result =
        std::make_shared<siconos::algebra::SimpleMatrix>(graph.vertices_number(), 25);

    int data_size = 0;
    for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend; ++vi) {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      /* create a visitor for specified classes */
      using ContactPointInspector = siconos::internal::Visitor<
          siconos::internal::Classes<
              siconos::modeling::NewtonEuler1DR, siconos::modeling::NewtonEuler3DR,
              siconos::modeling::NewtonEuler5DR, siconos::modeling::Lagrangian2d2DR,
              siconos::modeling::Lagrangian2d3DR,
              siconos::collision::native::bodies::CircleCircleR,
              siconos::collision::native::bodies::DiskDiskR,
              siconos::collision::native::bodies::DiskPlanR>,
          ContactPointVisitor>::Make;

      auto inspector = std::make_shared<ContactPointInspector>();
      inspector->inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      siconos::algebra::SiconosVector& data = inspector->answer;
      data_size = data.size();

      if (data_size == 0) {
        // Nothing is done since the relation does not appear as a relation
        // related to a contact points (perhaps a joint)
      } else {
        // We add at the end the number of ds1 and ds2
        data.resize(data_size + 2);
        DEBUG_EXPR(data.display(););
        auto& ds1 = *graph.properties(*vi).source;
        auto& ds2 = *graph.properties(*vi).target;

        data.setValue(data_size, ds1.number());
        data.setValue(data_size + 1, ds2.number());
        DEBUG_EXPR(data.display(););
        if (result->size(1) != data.size()) {
          result->resize(graph.vertices_number(), data.size());
        }
        result->setRow(current_row++, data);
        data_size += 2;
      }
    }
    result->resize(current_row, data_size);
    DEBUG_EXPR(result->display(););
    return result;
  }

  return nullptr;
}

/* Get contact informations */
/* default: a visitor that do nothing */
struct ContactInfoVisitor : public siconos::internal::SiconosVisitor {
  std::shared_ptr<siconos::modeling::Interaction> inter{nullptr};
  // std::vector<int> answer; better with a vector of int
  siconos::algebra::SiconosVector answer;

  template <typename T>
  void operator()(const T& rel) {}
};

/* then specializations : */
template <>
void ContactInfoVisitor::operator()(const siconos::modeling::NewtonEuler3DR& rel) {
  auto id = inter->number();
  answer.resize(10);
  // answer[0]= id;
  // answer[1]= 0; // reserve for ds1.number
  // answer[2]= 0; // reserve for ds2.number
  answer.setValue(0, id);
  answer.setValue(1, 0);
  answer.setValue(2, 0);
}

template <>
void ContactInfoVisitor::operator()(const siconos::collision::ContactR& rel) {
  auto id = inter->number();
  answer.resize(4);
  answer.setValue(0, id);
  answer.setValue(1, 0);
  answer.setValue(2, 0);
  if (rel.bodyShapeRecordB->staticBody) {
    answer.setValue(3, rel.bodyShapeRecordB->staticBody->number);
  } else
    answer.setValue(3, 0);
}

template <>
void ContactInfoVisitor::operator()(const siconos::collision::Contact5DR& rel) {
  auto id = inter->number();
  answer.resize(4);
  answer.setValue(0, id);
  answer.setValue(1, 0);
  answer.setValue(2, 0);
  if (rel.bodyShapeRecordB->staticBody) {
    answer.setValue(3, rel.bodyShapeRecordB->staticBody->number);
  } else
    answer.setValue(3, 0);
}

template <>
void ContactInfoVisitor::operator()(const siconos::collision::Contact2dR& rel) {
  auto id = inter->number();
  answer.resize(4);
  answer.setValue(0, id);
  answer.setValue(1, 0);
  answer.setValue(2, 0);
  if (rel.bodyShapeRecordB->staticBody) {
    answer.setValue(3, rel.bodyShapeRecordB->staticBody->number);
  } else
    answer.setValue(3, 0);
}

template <>
void ContactInfoVisitor::operator()(const siconos::collision::Contact2d3DR& rel) {
  auto id = inter->number();
  answer.resize(4);
  answer.setValue(0, id);
  answer.setValue(1, 0);
  answer.setValue(2, 0);
  if (rel.bodyShapeRecordB->staticBody) {
    answer.setValue(3, rel.bodyShapeRecordB->staticBody->number);
  } else
    answer.setValue(3, 0);
}

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::io::MechanicsIO::contactInfo(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds, unsigned int index_set) const {
  DEBUG_BEGIN("MechanicsIO::contactInfo");

  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  if (nsds.topology()->numberOfIndexSet() > 0) {
    auto& graph = *nsds.topology()->indexSet(index_set);
    unsigned int current_row;
    auto result = std::make_shared<siconos::algebra::SimpleMatrix>(graph.vertices_number(), 4);

    int data_size = 0;
    for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend; ++vi) {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      /* create a visitor for specified classes */
      using ContactInfoInspector = siconos::internal::Visitor<
          siconos::internal::Classes<
              siconos::modeling::NewtonEuler3DR, siconos::collision::ContactR,
              siconos::collision::Contact5DR, siconos::collision::Contact2dR,
              siconos::collision::Contact2d3DR>,
          ContactInfoVisitor>::Make;

      auto inspector = std::make_shared<ContactInfoInspector>();
      inspector->inter = graph.bundle(*vi);
      graph.bundle(*vi)->relation()->accept(inspector);
      auto& data = inspector->answer;
      data_size = data.size();

      if (data_size == 0) {
        // Nothing is done since the relation does not appear as a relation
        // related to a contact points (perhaps a joint)
      } else {
        // We add at the end the number of ds1 and ds2
        DEBUG_EXPR(data.display(););
        auto& ds1 = *graph.properties(*vi).source;
        auto& ds2 = *graph.properties(*vi).target;
        data.setValue(1, ds1.number());
        data.setValue(2, ds2.number());
      }
      if (result->size(1) != data.size()) {
        result->resize(graph.vertices_number(), data.size());
      }
      result->setRow(current_row++, data);
    }
    result->resize(current_row, data_size);
    DEBUG_EXPR(result->display(););
    return result;
  }
  DEBUG_END("SP::SimpleMatrix MechanicsIO::contactInfo");

  return nullptr;
}

namespace siconos::io {
/* Get contact work information */
/* default: a visitor that do nothing */

struct ContactContactWorkVisitor : public siconos::internal::SiconosVisitor {
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
  auto id = inter->number();
  answer.resize(6);
}

static void compute_contact_work_and_status(
    std::shared_ptr<siconos::modeling::Interaction> inter, double omega, double tol,
    siconos::algebra::SiconosVector& answer) {
  auto mu = siconos::internal::ask<ForMu>(*inter->nonSmoothLaw());
  auto e = siconos::internal::ask<ForE>(*inter->nonSmoothLaw());
  // Compute normal contact work
  auto vn_minus = inter->y_k(1).getValue(0);
  auto vn_plus = inter->y(1)->getValue(0);
  auto pn = inter->lambda(1)->getValue(0);

  double vn_average = omega * vn_plus + (1. - omega) * vn_minus;
  auto normal_contact_work = 0.5 * (vn_minus + vn_plus) * pn;
  answer.setValue(1, normal_contact_work);

  // Compute tangent contact work of impulse

  auto vt_1_minus = inter->y_k(1).getValue(1);
  auto vt_2_minus = inter->y_k(1).getValue(2);
  auto vt_1_plus = inter->y(1)->getValue(1);
  auto vt_2_plus = inter->y(1)->getValue(2);

  double vt_1_average = omega * vt_1_plus + (1. - omega) * vt_1_minus;
  double vt_2_average = omega * vt_2_plus + (1. - omega) * vt_2_minus;

  double pt_1 = inter->lambda(1)->getValue(1);
  double pt_2 = inter->lambda(1)->getValue(2);

  double tangent_contact_work = vt_1_average * pt_1 + vt_2_average * pt_2;
  answer.setValue(2, tangent_contact_work);

  // Compute directly work dissipated by friction impulse
  double norm_vt_average = sqrt(vt_1_average * vt_1_average + vt_2_average * vt_2_average);
  double friction_dissipation = mu * norm_vt_average * pn;
  answer.setValue(3, friction_dissipation);

  /* Compute contact status
   * Warning the status are consistent for the sticking and sliding
   * only with fully implicit discretization o NewtonImpact law
   * and not wih Fremond impact law
   */

  double norm_pt = sqrt(pt_1 * pt_1 + pt_2 * pt_2);
  double norm_vt_plus = sqrt(vt_1_plus * vt_1_plus + vt_2_plus * vt_2_plus);
  double norm_vt_minus = sqrt(vt_1_minus * vt_1_minus + vt_2_minus * vt_2_minus);
  if ((pn < tol) and (vn_plus + e * vn_minus > tol))
    answer.setValue(4, 0);  // take-off = 0
  else if ((pn > tol) and (vn_plus + e * vn_minus < tol)) {
    if ((norm_pt - mu * pn > tol)) {
      // std::cout << "WARNING: the impulse is outside the Coulomb cone  " << std::endl;
      answer.setValue(4, -3);  // outside the cone = -3
    } else if ((norm_pt - mu * pn < -tol)) {
      // std::cout << "the impulse is in the *interior* of  the Coulomb cone  " << std::endl;
      // std::cout << "norm_vt_plus  " << norm_vt_plus << std::endl;
      if (norm_vt_plus > tol) {
        // std::cout << "WARNING: but the norm of vt is not zero  " << std::endl;
        answer.setValue(4, -2);  // sticking with a non zero slifing velocity = -2
      }
      answer.setValue(4, 1);  // sticking = 1
    } else {
      // std::cout << "the impulse is on the *boundary* of the Coulomb cone  " << std::endl;
      // std::cout << "norm_vt_plus  " << norm_vt_plus << std::endl;
      answer.setValue(4, 2);  // sliding = 2
    }
  } else
    answer.setValue(4, -1);  // undetermined = -1

  if ((pn > tol) and (vn_minus > tol)) {
    // std::cout << "WARNING: we apply the impact law of positive velocity " << std::endl;
    // std::cout << "pn " << pn << " vn minus " << vn_minus << " vn plus " << vn_plus
    // 		<< " normal_contact_work " << normal_contact_work
    // 		<< " -e * vn_minus   " << -e*vn_minus
    // 		<< std::endl;
    answer.setValue(5, normal_contact_work);
  }
  // double id = inter->number();
  // std::cout << "\nid "<< id << std::endl;
  // std::cout << " e "<< e  << " mu "<< mu << std::endl;
  // std::cout << " tol "<< tol << " omega " << omega << std::endl;
  // std::cout << "vn_plus "<< vn_plus << std::endl;
  // std::cout << "vn_minus "<< vn_minus << std::endl;
  // std::cout << "pn "<< pn << std::endl;
  // std::cout << "normal_contact_work  "<< normal_contact_work  << std::endl;

  // std::cout << "vt_plus "<< vt_1_plus << " " << vt_2_plus <<  std::endl;
  // std::cout << "vt_minus "<< vt_1_minus << " " << vt_2_minus <<  std::endl;
  // std::cout << "pt "<< pt_1 << " "  << pt_2 << std::endl;
  // std::cout << "tangent_contact_work  "<< tangent_contact_work  << std::endl;

  // std::cout << "friction_dissipation  "<< friction_dissipation << std::endl;

  // std::cout << "norm_pt  "<< norm_pt  << std::endl;
  // std::cout << "norm_pt - mu* pn  "<< norm_pt -mu *pn   << std::endl;
  // std::cout << "vn_plus + e * vn_minus  " << vn_plus + e * vn_minus   << std::endl;
  // std::cout << "status   "<<   answer.getValue(4) << std::endl;
}
}  // namespace siconos::io

template <>
void siconos::io::ContactContactWorkVisitor::operator()(
    const siconos::collision::ContactR& rel) {
  auto id = inter->number();
  answer.resize(6);
  answer.setValue(0, id);
  compute_contact_work_and_status(inter, omega, tol, answer);
}

template <>
void siconos::io::ContactContactWorkVisitor::operator()(
    const siconos::collision::Contact5DR& rel) {
  auto id = inter->number();
  answer.resize(6);
  answer.setValue(0, id);

  compute_contact_work_and_status(inter, omega, tol, answer);
}

template <>
void siconos::io::ContactContactWorkVisitor::operator()(
    const siconos::collision::Contact2dR& rel) {
  auto id = inter->number();
  answer.resize(6);
  answer.setValue(0, id);

  compute_contact_work_and_status(inter, omega, tol, answer);
}

template <>
void siconos::io::ContactContactWorkVisitor::operator()(
    const siconos::collision::Contact2d3DR& rel) {
  auto id = inter->number();
  answer.resize(6);
  answer.setValue(0, id);

  compute_contact_work_and_status(inter, omega, tol, answer);
}

std::shared_ptr<siconos::algebra::SimpleMatrix> siconos::io::MechanicsIO::contactContactWork(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds, unsigned int index_set,
    double omega, double tol) const {
  DEBUG_BEGIN("SimpleMatrix MechanicsIO::contactContactWork");

  siconos::graphs::InteractionsGraph::VIterator vi, viend;
  if (nsds.topology()->numberOfIndexSet() > 0) {
    auto& graph = *nsds.topology()->indexSet(index_set);

    unsigned int current_row;
    auto result =
        std::make_shared<siconos::algebra::SimpleMatrix>(graph.vertices_number(), 25);

    int data_size = 0;
    for (current_row = 0, std::tie(vi, viend) = graph.vertices(); vi != viend; ++vi) {
      DEBUG_PRINTF("process interaction : %p\n", &*graph.bundle(*vi));

      /* create a visitor for specified classes */
      using ContactContactWorkInspector = siconos::internal::Visitor<
          siconos::internal::Classes<
              siconos::modeling::NewtonEuler3DR, siconos::collision::ContactR,
              siconos::collision::Contact5DR, siconos::collision::Contact2dR,
              siconos::collision::Contact2d3DR>,
          ContactContactWorkVisitor>::Make;
      auto inspector = std::make_shared<ContactContactWorkInspector>();
      inspector->inter = graph.bundle(*vi);
      inspector->tol = tol;
      inspector->omega = omega;
      graph.bundle(*vi)->relation()->accept(inspector);
      auto& data = inspector->answer;
      data_size = data.size();

      if (data_size == 0) {
        // Nothing is done since the relation does not appear as a relation
        // related to a contact points (perhaps a joint)
      } else {
      }
      if (result->size(1) != data.size()) {
        result->resize(graph.vertices_number(), data.size());
      }
      result->setRow(current_row++, data);
    }
    result->resize(current_row, data_size);
    DEBUG_EXPR(result->display(););
    return result;
  }
  DEBUG_END("SP::SimpleMatrix MechanicsIO::contactContactWork");

  // result->display();
  return nullptr;
}
