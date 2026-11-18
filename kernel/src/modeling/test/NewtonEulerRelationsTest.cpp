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
#include "NewtonEulerRelationsTest.hpp"

#include <memory>
#include <random>
#include <sstream>

#include "RotationQuaternion.hpp"
#include "SiconosKernel.hpp"

using siconos::algebra::BlockVector;
using siconos::algebra::SiconosMatrix;
using siconos::algebra::SiconosMatrix33;
using siconos::algebra::SiconosVector;
using siconos::algebra::SiconosVector3;
using siconos::algebra::SiconosVector6;
using siconos::algebra::SiconosVector7;
using namespace siconos::modeling;

CPPUNIT_TEST_SUITE_REGISTRATION(NewtonEulerRelationsTest);

void NewtonEulerRelationsTest::setUp() {}

void NewtonEulerRelationsTest::tearDown() {}

namespace {

constexpr int nTrials = 20;
constexpr double tol = 1e-12;

template <class Vec>
Vec randomVector(std::mt19937& gen) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  Vec v;
  for (Eigen::Index i = 0; i < v.size(); ++i) v(i) = dist(gen);
  return v;
}

SiconosVector7 randomPose(std::mt19937& gen) {
  SiconosVector7 q = randomVector<SiconosVector7>(gen);
  q.tail<4>().normalize();
  return q;
}

/** Builds the relation and the interaction for one or two random bodies, computes
 *  H_NE_prod_T and returns || H_NE_prod_T * twists - y_ref ||.
 */
template <class Rel>
double jacobianError(std::shared_ptr<NonSmoothLaw> nslaw, bool twoBodies, std::mt19937& gen) {
  const auto ySize = nslaw->size();
  auto rel = std::make_shared<Rel>();
  Interaction inter(nslaw, rel);
  inter.setHas2Bodies(twoBodies);
  inter.setDSSizes(twoBodies ? 12 : 6);

  const SiconosVector7 q1 = randomPose(gen);
  const SiconosVector7 q2 = randomPose(gen);
  const SiconosMatrix33 inertia = SiconosMatrix33::Identity();
  NewtonEulerDS ds1(q1, SiconosVector6::Zero(), 1.0, inertia, siconos::algebra::copy_t);
  NewtonEulerDS ds2(q2, SiconosVector6::Zero(), 1.0, inertia, siconos::algebra::copy_t);
  inter.initialize_read_access_to_ds_variables(ds1, twoBodies ? ds2 : ds1);

  // Contact points (in the body frames) and normal (in the frame of body 2, absolute if
  // there is only one body)
  rel->setRelPc1(randomVector<SiconosVector3>(gen));
  rel->setRelPc2(randomVector<SiconosVector3>(gen));
  rel->setRelNc(randomVector<SiconosVector3>(gen).normalized());
  SiconosVector y(ySize);
  if (twoBodies) {
    const SiconosVector q2dyn = q2;
    rel->computehFromRelativeContactPoints(0., q1, q2dyn, y);
  } else {
    rel->computehFromRelativeContactPoints(0., q1, std::nullopt, y);
  }

  BlockVector q0;
  q0.insertPtr(std::make_shared<SiconosVector>(q1));
  if (twoBodies) q0.insertPtr(std::make_shared<SiconosVector>(q2));
  rel->computeH_NE_prod_T(inter, q0);
  const SiconosMatrix H = rel->H_NE_prod_T();

  // Twists: translational velocity (absolute frame), angular velocity (body frame)
  const SiconosVector6 twist1 = randomVector<SiconosVector6>(gen);
  const SiconosVector6 twist2 = randomVector<SiconosVector6>(gen);
  SiconosVector twists(twoBodies ? 12 : 6);
  twists.head<6>() = twist1;
  if (twoBodies) twists.tail<6>() = twist2;
  const SiconosVector yH = H * twists;

  // Reference: relative velocity at pc1 in the contact frame
  SiconosMatrix33 R1, R2;
  siconos::geometry::computeRotationMatrix(q1, R1);
  siconos::geometry::computeRotationMatrix(q2, R2);
  const SiconosVector3 pc = rel->pc1();
  SiconosVector3 n = rel->nc(), t1, t2;
  CPPUNIT_ASSERT(siconos::geometry::orthoBaseFromVector(n, t1, t2));

  const SiconosVector3 omega1 = R1 * twist1.tail<3>();
  const SiconosVector3 omega2 = R2 * twist2.tail<3>();
  SiconosVector3 u = twist1.head<3>() + omega1.cross(pc - q1.head<3>());
  SiconosVector3 omegaRel = omega1;
  if (twoBodies) {
    u -= twist2.head<3>() + omega2.cross(pc - q2.head<3>());
    omegaRel -= omega2;
  }

  SiconosVector yRef(ySize);
  yRef(0) = n.dot(u);
  if (ySize >= 3) {
    yRef(1) = t1.dot(u);
    yRef(2) = t2.dot(u);
  }
  if (ySize == 5) {
    yRef(3) = t1.dot(omegaRel);
    yRef(4) = t2.dot(omegaRel);
  }
  return (yH - yRef).norm();
}

template <class Rel>
void checkJacobian(const std::string& name, std::shared_ptr<NonSmoothLaw> nslaw) {
  std::mt19937 gen(42);
  for (bool twoBodies : {false, true}) {
    for (int trial = 0; trial < nTrials; ++trial) {
      const double err = jacobianError<Rel>(nslaw, twoBodies, gen);
      std::ostringstream msg;
      msg << name << " (" << (twoBodies ? "two bodies" : "one body") << ", trial " << trial
          << "): || H_NE_prod_T * twist - u_ref || = " << err;
      CPPUNIT_ASSERT_MESSAGE(msg.str(), err < tol);
    }
  }
}

}  // namespace

void NewtonEulerRelationsTest::testNewtonEuler1DRJacobian() {
  checkJacobian<NewtonEuler1DR>("NewtonEuler1DR", std::make_shared<NewtonImpactNSL>(0.));
}

void NewtonEulerRelationsTest::testNewtonEuler3DRJacobian() {
  checkJacobian<NewtonEuler3DR>("NewtonEuler3DR",
                                std::make_shared<NewtonImpactFrictionNSL>(0., 0., 0.5, 3));
}

void NewtonEulerRelationsTest::testNewtonEulerVelocityAngularVelocityRJacobian() {
  checkJacobian<NewtonEulerVelocityAngularVelocityR>(
      "NewtonEulerVelocityAngularVelocityR",
      std::make_shared<NewtonImpactRollingFrictionNSL>(0., 0., 0.5, 0.1, 5));
}
