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
#ifndef HARMONICBC_HPP
#define HARMONICBC_HPP

#include "BoundaryCondition.hpp"
#include "SiconosVector.hpp"

namespace siconos::modeling {
/** \brief This class models a simple harmonic boundary conditions for
 *   prescribing the velocities in a Dynamical System.
 *
 * A simple boundary condition is considered to fix a component \f$ j \f$ of
 *   the velocity vector, i.e.,
 *  \f$ v_j(t) = a +  b cos( \omega t+ \phi) \f$.
 *
 */
class HarmonicBC : public BoundaryCondition {
 public:
  /** Constructor
   * \param newVelocityIndices the indices of the velocity subjected to prescribed velocities
   * \param a constant value for additive term of the prescribed velocity
   * \param b constant value for multiplicative term of the prescribed velocity
   * \param omega frequency
   * \param phi phase
   */
  HarmonicBC(Indices&& newVelocityIndices, double a, double b, double omega, double phi)
      : BoundaryCondition(std::move(newVelocityIndices)),
        aCoeff_(a),
        bCoeff_(b),
        omega_(omega),
        phi_(phi) {};

  HarmonicBC(const Indices& newVelocityIndices,
             Eigen::Ref<siconos::algebra::SiconosVector> newa,
             Eigen::Ref<siconos::algebra::SiconosVector> newb,
             Eigen::Ref<siconos::algebra::SiconosVector> omega,
             Eigen::Ref<siconos::algebra::SiconosVector> phi);

  /** destructor */
  virtual ~HarmonicBC() noexcept = default;

  /** default function to compute the precribed velocities
   *
   *  \param  time : the current time
   */
  virtual void computePrescribedVelocity(double time);

 protected:
  ACCEPT_SERIALIZATION(HarmonicBC);

  /** Constant additive term of the prescribed velocity  */
  double aCoeff_ = 0.;
  /** Constant multiplicative term of the prescribed velocity  */
  double bCoeff_ = 0;
  ;
  /** Constant frequency  */
  double omega_ = 0;
  ;
  /** Constant phase  */
  double phi_ = 0;
  ;

  /** Constant additive term of the prescribed velocity  */
  std::unique_ptr<siconos::algebra::MapVectorType> a_view_{nullptr};
  /** Constant multiplicative term of the prescribed velocity  */
  std::unique_ptr<siconos::algebra::MapVectorType> b_view_{nullptr};
  /** Constant frequency  */
  std::unique_ptr<siconos::algebra::MapVectorType> omega_view_{nullptr};
  /** Constant phase  */
  std::unique_ptr<siconos::algebra::MapVectorType> phi_view_{nullptr};
};
}  // namespace siconos::modeling
#endif  // HARMONICBC_HPP
