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

#include "TransportCableProfile.h"

#include <random>

#include "CableTools.h"
#include "Carriers.h"
#include "MechanicalProperties.h"
#include "PulleyWrapping.h"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "TransportCableModel.h"
#include "TransportCableResult.h"

siconos::fem::cable::TransportCableProfile::TransportCableProfile(
    const TransportCableModel &a_model, TransportCableResult &a_results)
    : r_model(a_model), r_results(a_results) {}

void siconos::fem::cable::TransportCableProfile::computeInitialProfile(int nb_nodes,
                                                                       double a_tol,
                                                                       int a_nmax) {
  // calcul les positions, tensions des cordes

  // Get mechanical properties and carriers from the model
  MechanicalProperties meca = r_model.mechanicalProperties();
  const auto &vehicles = r_model.get_carriers();
  // Update rho as the sum of both cable and vehicles rho
  meca.set_rho(meca.linearDensity() + vehicles.linearDensity);

  // Build all ropes and compute their initial profiles (Catenary)
  // up ...
  r_results.ropes_up.computeCatenary(meca, r_model.get_pylons_up(), nb_nodes, a_tol, a_nmax);
  meca.set_T(r_results.ropes_up.initialTension());

  // and down
  r_results.ropes_down.computeCatenary(meca, r_model.get_pylons_down(), nb_nodes, a_tol,
                                       a_nmax);

  // Build supports
  r_results.prepareSupport();
}

void siconos::fem::cable::TransportCableProfile::initializeFEM(int nb_elem, double a_eps,
                                                               double a_tol) {
  // -- Computes the lengths of the cable around the pulleys (stations)
  assert(r_results.topPulleyId >= 0);
  assert(r_results.downPulleyId >= 0);
  auto topPulley = std::dynamic_pointer_cast<PulleyWrapping>(
      r_results.supports[static_cast<size_t>(r_results.topPulleyId)]);

  auto downPulley = std::dynamic_pointer_cast<PulleyWrapping>(
      r_results.supports[static_cast<size_t>(r_results.downPulleyId)]);

  double topLength = topPulley->length(r_results.ropes_down);
  double downLength = downPulley->length(r_results.ropes_up);

  // -- Computes the total length of the cable --
  //  (pulleys'wrappers + top and down ropeways lengths)
  r_results.totalLength =
      r_results.ropes_up.length() + topLength + r_results.ropes_down.length() + downLength;

  // -- Computes the number of elements in each part --
  // In the wrappings around the pulleys:
  int nElemInTopPulley = static_cast<int>(rint(nb_elem * topLength / r_results.totalLength));
  int nElemInDownPulley = static_cast<int>(rint(nb_elem * downLength / r_results.totalLength));

  double elementLength = r_results.totalLength / nb_elem;
  int nElemUpCable =
      r_results.ropes_up.computeNumberOfElements(elementLength, r_results.totalLength);
  int nElemDownCable =
      r_results.ropes_down.computeNumberOfElements(elementLength, r_results.totalLength);
  r_results.numberOfElements =
      nElemInTopPulley + nElemUpCable + nElemInDownPulley + nElemDownCable;

  // -- Initialize positions, tension and reaction vectors --
  auto &q = r_results.q;
  q.resize(3 * r_results.numberOfElements);
  auto &R = r_results.R;
  R.resize(3 * r_results.numberOfElements);
  auto &TS = r_results.TS;
  TS.resize(r_results.numberOfElements);

  // Fills q, R and TS corresponding to the ropes in the up cable.
  // --> solves catenary for the FE mesh
  // Warning: last node is not saved in q
  // The offset is used to update the starting position in global vector q, R and TS
  // where results must be computed.
  // up rope
  auto offset = r_results.ropes_up.initializeFEM(q, R, TS, 0);
  // top pullet
  offset = topPulley->computeMesh(nElemInTopPulley + 1, q, offset);
  // down rope
  offset = r_results.ropes_down.initializeFEM(q, R, TS, offset);
  // down pulley
  offset = downPulley->computeMesh(nElemInDownPulley + 1, q, offset);
  assert(offset == r_results.numberOfElements);

  // Computes the real elements size
  r_results.elementLength = r_results.totalLength / r_results.numberOfElements;

  // Applies weights due to carriers
  r_results.weightVector = distribute_carriers_weight(
      r_model.get_carriers(), r_results.numberOfElements, r_results.totalLength);

  //  compute_ineq_constraint(q, a_tol);
  r_results.prepareIneqConstraint(q.size() / 3);
  // computeConstraintsSparse(q, a_tol, r_results.gVector, r_results.jacobian_g_Over_q);
  computeConstraints(q, a_tol, r_results.gVector, r_results.jacobian_g_Over_q);
  /*dgi = np.zeros(nb_node)
  dgi[gi <= 0] = gi[gi <= 0]
  q = q - np.matmul(np.transpose(Gi), 1.1*dgi)*/
  // const auto &gVector = r_results.gVector;
  double k = 1. + a_eps;
  q -= k * r_results.jacobian_g_Over_q * r_results.gVector;

  // for (auto i = 0; i < r_results.numberOfElements; i++) {
  //   siconos::algebra::SiconosVector3 point;
  //   point.setZero();
  //   if (gVector(0) < 0) {
  //     for (auto j = 0; j < r_results.numberOfElements; j++) {
  //       point += jacobian->row(j).segment<3>(i);
  //     }
  //     point *= k * gVector(i);
  //     q.segment<3>(i) -= point;
  //   }
  // }
}

void siconos::fem::cable::TransportCableProfile::compute_ineq_constraint(
    const siconos::algebra::SiconosVector &a_X, double a_tol) {
  /*
  @author: charl

  get_ineq_constraint(q,supports)

  where

  X is the coordinates of cable particles
  supports is the list of obstacles which i-th element is containing [ positions, radius ]
  associated to piles (cylinder with circled base in xy plane) pulleys is the list of pulleys
  which i-th element is containing [ positions, radius ] associated to pulleys (cylinder with
  circled base in xz plane)

  Return the inequality constraint vector for the support contained into supports and the
  active set vector which i-th element is :
          - NaN if the constraint is inactive
          - Obstacle index in supports if the constraint is active
  */
  // r_results.prepareIneqConstraint((int)a_X.size());

  // auto &c = r_results.contacts;
  // auto &gVector = r_results.gVector;
  // auto jacobian_g_Over_q = r_results.jacobian_g_Over_q;
  // auto T = r_results.T;
  // auto &supports = r_results.supports;

  // auto nb_nodes =
  //     a_X.size() / 3;  // TODO get this elsewhere ... This should be numberOfElements?

  // for (auto &s : supports) {
  //   for (int pos = 0; pos < nb_nodes; ++pos) {
  //     s->compute(a_X.segment<3>(pos), a_tol, gVector(pos),
  //                jacobian_g_Over_q.col(pos).segment<3>(pos), T->col(pos).segment<3>(pos),
  //                c[pos]);
  //   }
  // }
}

void siconos::fem::cable::TransportCableProfile::computeConstraints(
    const siconos::algebra::SiconosVector &cableNodesPositions, double tol,
    Eigen::Ref<siconos::algebra::SiconosVector> distances,
    Eigen::Ref<siconos::algebra::SiconosDenseMatrix> jacobian) {
  const int nb = cableNodesPositions.size() / 3;
  for (auto &s : r_results.supports) {
    const auto &center = s->center();
    auto radius = s->radius();
    if (auto pulley = std::dynamic_pointer_cast<PulleyWrapping>(s)) {
      for (int i = 0; i < nb; ++i) {
        double dx = cableNodesPositions(3 * i) - center(0);
        double dy = cableNodesPositions(3 * i + 1) - center(1);
        double norm = std::sqrt(dx * dx + dy * dy);

        if ((norm - radius) <= tol) {
          jacobian(3 * i, i) = dx / norm;
          jacobian(3 * i + 1, i) = dy / norm;
          distances(i) = norm - radius;
        }
      }
    } else {
      for (int i = 0; i < nb; ++i) {
        double dx = cableNodesPositions(3 * i) - center(0);
        double dz = cableNodesPositions(3 * i + 2) - center(2);
        double norm = std::sqrt(dx * dx + dz * dz);

        if ((norm - radius) <= tol) {
          jacobian(3 * i, i) = dx / norm;
          jacobian(3 * i + 2, i) = dz / norm;
          distances(i) = norm - radius;
        }
      }
    }
  }
  for (int i = 0; i < nb; ++i) {
    if (distances(i) > tol) distances(i) = 1.;
  }
}

void siconos::fem::cable::TransportCableProfile::computeConstraintsSparse(
    const siconos::algebra::SiconosVector &cableNodesPositions, double tol,
    Eigen::Ref<siconos::algebra::SiconosVector> distances,
    siconos::algebra::SiconosSparseMatrix &jacobian) {
  const siconos::algebra::Index nb = cableNodesPositions.size() / 3;
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(2 * siconos::algebra::to_unsigned<size_t>(nb));  // max 2 per col¨
  for (auto &s : r_results.supports) {
    const auto &center = s->center();
    auto radius = s->radius();
    if (auto pulley = std::dynamic_pointer_cast<PulleyWrapping>(s)) {
      for (siconos::algebra::Index i = 0; i < nb; ++i) {
        double dx = cableNodesPositions(3 * i) - center(0);
        double dy = cableNodesPositions(3 * i + 1) - center(1);
        double norm = std::sqrt(dx * dx + dy * dy);

        if ((norm - radius) <= tol) {
          triplets.emplace_back(3 * i, i, dx / norm);
          triplets.emplace_back(3 * i + 1, i, dy / norm);
          distances(i) = norm - radius;
        }
      }
    } else {
      for (siconos::algebra::Index i = 0; i < nb; ++i) {
        double dx = cableNodesPositions(3 * i) - center(0);
        double dz = cableNodesPositions(3 * i + 2) - center(2);
        double norm = std::sqrt(dx * dx + dz * dz);

        if ((norm - radius) <= tol) {
          triplets.emplace_back(3 * i, i, dx / norm);
          triplets.emplace_back(3 * i + 2, i, dz / norm);
          distances(i) = norm - radius;
        }
      }
    }
  }
  for (int i = 0; i < nb; ++i) {
    if (distances(i) > tol) distances(i) = 1.;
  }
  jacobian.setFromTriplets(triplets.begin(), triplets.end());
  // jacobian.makeCompressed();
}

// Free functions

siconos::algebra::SiconosVector siconos::fem::cable::distribute_carriers_weight(
    const Carriers &vehicles, int nb_elem, double totalLength) {
  siconos::algebra::SiconosVector weight;
  weight.resize(nb_elem);
  weight.setZero();

  double distanceBetweenCarriers = vehicles.distanceBetweenCarriers;
  assert(distanceBetweenCarriers > 0 &&
         "Distance between carriers must be strictly positive.");

  // Estimated number of carriers
  int nb_vehicles = static_cast<int>(totalLength / distanceBetweenCarriers);

  //  mass of a carrier
  double mass = vehicles.linearDensity * totalLength / nb_vehicles;

  // If not set properly, the first carrier pisition is
  //  set randomly between 0 and distanceBetweenCarriers on the cable
  double initialPosition = vehicles.firstCarrierInitialPosition;
  double start = initialPosition;
  if (initialPosition <= 0.) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, distanceBetweenCarriers);
    start = dis(gen);
  }

  // nb_elem points distributed between 0 and totalLength
  auto lc = siconos::fem::cable::tools::linspace(0., totalLength, nb_elem);
  lc.erase(lc.begin());
  // Utility vector to find the index corresponding
  // to each vehicle
  std::vector<int> ind(lc.size());
  for (size_t i = 0; i < lc.size(); ++i)
    ind[i] = static_cast<int>((lc[i] - start) / distanceBetweenCarriers);

  // Find carriers and set mass to each of them

  for (int k = 1; k < nb_vehicles; ++k) {
    // --> search in ind the first occurence of each k
    // if found, it corresponds to a vehicle (or is close to ...)
    auto it = std::find_if(ind.begin(), ind.end(), [k](int val) { return val > k; });
    if (it != ind.end()) {
      siconos::algebra::Index index =
          siconos::algebra::to_index(std::distance(ind.begin(), it));
      weight(index) = mass;
    }
  }
  return weight;  // RVO
}
