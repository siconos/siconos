// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file VonMisesPlasticityNSL.hpp
  Von Mises (J2) plasticity non-smooth law
*/

#ifndef VONMISESPLASTICITYNSL_H
#define VONMISESPLASTICITYNSL_H

#include "NonSmoothLaw.hpp"

namespace siconos::modeling {

/** Von Mises (J2) plasticity law with associated flow rule, for 2D stresses.
 *
 * The output y of the interaction is a stress in Mandel notation, y = (σxx, σyy, √2 σxy),
 * and must stay in the elastic domain K; the multiplier λ is (minus h times) a plastic
 * strain rate in the normal cone (associated flow rule): 0 ∈ λ + N_K(y), also in Mandel
 * notation (ε̇xx, ε̇yy, √2 ε̇xy).
 *
 * - plane stress (σzz = 0): K = {σxx² + σyy² − σxx σyy + 3 σxy² ≤ σ_y²};
 * - plane strain (incompressible plastic flow, σzz = (σxx + σyy)/2):
 *   K = {(σxx − σyy)²/4 + σxy² ≤ σ_y²/3}.
 *
 * Solved by fem::nonsmooth_formulations::Plasticity with the numerics PlasticityProblem
 * (Von Mises model, Mandel coordinates). See docs/devel/FEM_PLASTICITY_SCHEME.md.
 */
class VonMisesPlasticityNSL : public NonSmoothLaw {
 public:
  /** 2D stress hypothesis */
  enum class Hypothesis { plane_strain, plane_stress };

 private:
  ACCEPT_SERIALIZATION(VonMisesPlasticityNSL);

  /** yield stress */
  double _sigma_y{0.};

  /** 2D hypothesis */
  Hypothesis _hypothesis{Hypothesis::plane_stress};

  /** default constructor */
  VonMisesPlasticityNSL() = default;

 public:
  /** constructor
   *
   *  \param sigma_y yield stress (> 0)
   *  \param hypothesis plane stress (default) or plane strain
   *  \param size size of the law (3: Mandel components of a 2D stress)
   */
  VonMisesPlasticityNSL(double sigma_y, Hypothesis hypothesis = Hypothesis::plane_stress,
                        unsigned int size = 3)
      : NonSmoothLaw(size), _sigma_y(sigma_y), _hypothesis(hypothesis) {};

  /** Destructor */
  ~VonMisesPlasticityNSL() noexcept = default;

  /** not implemented (throws) */
  bool isVerified(void) const override;

  /** \return the yield stress */
  inline double sigma_y() const { return _sigma_y; };

  /** set the yield stress
   *  \param newVal the yield stress
   */
  inline void setSigma_y(double newVal) { _sigma_y = newVal; };

  /** \return the 2D hypothesis */
  inline Hypothesis hypothesis() const { return _hypothesis; };

  /** print the data to the screen */
  void display() const override;

  // visitors hook
  virtual void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // VONMISESPLASTICITYNSL_H
