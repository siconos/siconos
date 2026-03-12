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
/*! \file TransportCableProfile.h

*/

#pragma once

#include <nlohmann/json.hpp>

#include "Carriers.h"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::fem::cable {

class TransportCableModel;
class TransportCableResult;
class Point;

/**  */
class TransportCableProfile {
 public:
  TransportCableProfile(const TransportCableModel &a_model, TransportCableResult &a_results);

  virtual ~TransportCableProfile() noexcept = default;

  /** Build all ropes for each cable (up and down) and initialize (Catenary) their profiles
     + build and initialize supports

     \param[in] nb_nodes number of nodes by segment
     \param[in] tol tolerance used to compute initial profile
     \param[in] nmax max number of iterations used to compute initial profile

  */
  void computeInitialProfile(int nb_nodes, double tol = 1e-20, int nmax = 20);

  void initializeFEM(int nb_elem, double a_eps = 0.1, double a_tol = 1e-3);

 private:
  const TransportCableModel &r_model;
  TransportCableResult &r_results;

  /**
     \param a_X vector of positions, coordinates of cable 'particles'
     \param a_tol tolerance used to activate contacts
   */
  void compute_ineq_constraint(const siconos::algebra::SiconosVector &a_X,
                               double a_tol = 1e-3);

  void computeConstraints(const siconos::algebra::SiconosVector &cableNodesPositions,
                          double tol, Eigen::Ref<siconos::algebra::SiconosVector> distances,
                          Eigen::Ref<siconos::algebra::SiconosDenseMatrix> jacobian);
  void computeConstraintsSparse(const siconos::algebra::SiconosVector &cableNodesPositions,
                                double tol,
                                Eigen::Ref<siconos::algebra::SiconosVector> distances,
                                siconos::algebra::SiconosSparseMatrix &jacobian);
};
// Free functions

/** Compute the weight vector, due to mass distribution of carriers

    \param carriers the set of vehicles
    \param nb_elem number of elements in the mesh (i.e. size of weight vector)
    \param totalLength total length of the cable (up, down and stations)
    \return weight vector, with positive values at indices where vehicles are supposed to be
  */
siconos::algebra::SiconosVector distribute_carriers_weight(const Carriers &carriers,
                                                           int nb_elem, double totalLength);

}  // namespace siconos::fem::cable
