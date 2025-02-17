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

#ifndef MechanicsIO_hpp
#define MechanicsIO_hpp

#include <memory>
#include <optional>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

// #include <MechanicsFwd.hpp>

// #ifdef HAVE_SICONOS_MECHANISMS
// #include <MechanismsFwd.hpp>
// #endif
// #include <SiconosFwd.hpp>
// #include <SiconosPointers.hpp>

namespace siconos::modeling {

class NonSmoothDynamicalSystem;
}

namespace siconos::io {

struct GetPosition;
struct GetVelocity;
struct ForMu;
struct ForE;
struct ContactPointVisitor;
struct ContactPointDomainVisitor;
struct ContactContactWorkVisitor;

class MechanicsIO {
 protected:
  template <typename T, typename G>
  siconos::algebra::SiconosMatrix visitAllVerticesForVector(const G& graph) const;

  template <typename T, typename G>
  siconos::algebra::SiconosVector visitAllVerticesForDouble(const G& graph) const;

  MechanicsIO(const MechanicsIO&) = delete;
  MechanicsIO& operator=(const MechanicsIO&) = delete;
  MechanicsIO(MechanicsIO&&) = delete;
  MechanicsIO& operator=(MechanicsIO&&) = delete;

 public:
  /** default constructor
   */
  MechanicsIO() = default;
  ~MechanicsIO() noexcept = default;

  /** Collect  positions from all dynamical systems and save them in a matrix
   *
   *  \param nsds the nonsmooth dynamical system
   *  \return a matrix whith column(i) = [ds(i).number, ds(i).q_read()].T
   */
  siconos::algebra::SiconosMatrix positions(
      const siconos::modeling::NonSmoothDynamicalSystem& nsds) const;

  /** Collect velocities/twists from all dynamical systems and save them in a matrix
   *  \param nsds current nonsmooth dynamical system
   *  \return  a matrix whith column(i) = [ds(i).number, ds(i).velocity_read()].T
   */
  siconos::algebra::SiconosMatrix velocities(
      const siconos::modeling::NonSmoothDynamicalSystem& nsds) const;

  /** get the coordinates of all contact points, normals, reactions and velocities
   *  \param nsds current nonsmooth dynamical system
   *  \param index_set the index set number.
   *  \return a matrix where the columns are mu x y z, nx, ny, nz, rx, ry, rz, vx, vy, vz, ox,
   *   oy, oz, id
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> contactPoints(
      const siconos::modeling::NonSmoothDynamicalSystem& nsds,
      unsigned int index_set = 1) const;

  /** get the contact information that is the ds linked by the interaction
   *  \param nsds current nonsmooth dynamical system
   *  \param index_set the index set number.
   *  \return a matrix where the columns are interaction id, ds1 number, ds2 number, static
   *   object number (if possible)
   */

  std::optional<siconos::algebra::SiconosMatrix> contactInfo(
      const siconos::modeling::NonSmoothDynamicalSystem& nsds,
      unsigned int index_set = 1) const;

  /** get the dissipation values  of all contact points

      \param nsds current nonsmooth dynamical system
      \param index_set the index set number.
      \param omega the value of the weigth for the weight in the computaion of the contact work
      by default omega =1/2 and the contact work corresponds to the theoretical formula
      1/2 (v^+ + v^-)^\top p
      otherwise it corresponds to v_{k+omega} p
      \param tol double for the computation of contact status
      \return a matrix where the columns are id, normal contact work, tangent contact work,
      friction dissipation, contact status
  */
  std::shared_ptr<siconos::algebra::SiconosMatrix> contactContactWork(
      const siconos::modeling::NonSmoothDynamicalSystem& nsds, unsigned int index_set = 1,
      double omega = 0.5, double tol = 1e-08) const;

  /** get the domain of each contact point
   *  \param nsds current nonsmooth dynamical system
   *  \return a matrix where the columns are domain, id
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> domains(
      const siconos::modeling::NonSmoothDynamicalSystem& nsds) const;
};

}  // namespace siconos::io
#endif
