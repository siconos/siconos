// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file DruckerPragerPlasticityNSL.hpp
  Drucker-Prager (in-plane Mohr-Coulomb) plasticity non-smooth law
*/

#ifndef DRUCKERPRAGERPLASTICITYNSL_H
#define DRUCKERPRAGERPLASTICITYNSL_H

#include "NonSmoothLaw.hpp"

namespace siconos::modeling {

/** Drucker-Prager plasticity law for 2D stresses, with cohesion and dilatancy.
 *
 * The output y of the interaction is a stress in Voigt notation, y = (σxx, σyy, σxy)
 * (tension positive). With σ_N = -(σxx + σyy)/√2 and σ_T = ((σxx - σyy)/√2, √2 σxy)
 * (the cone coordinates of the numerics PlasticityProblem), the elastic domain is the cone
 *
 *   K = { ‖σ_T‖ ≤ η σ_N + c },
 *
 * i.e. R ≤ -η p + c/√2 with p = (σxx + σyy)/2 and R = √(((σxx - σyy)/2)² + σxy²) the
 * radius of the Mohr circle: with η = sin φ and c = √2 c₀ cos φ this is the Mohr-Coulomb
 * criterion of the plane (exact in plane strain when σzz is the intermediate principal
 * stress), see mohrCoulombPlaneStrain.
 *
 * Flow rule: -(w + θ ‖w_T‖ e_N) ∈ N_K(y), w = -h ε̇_p in cone coordinates (the multiplier
 * λ of the interaction): θ = 0 is the associated flow; θ = sin φ - sin ψ for a dilatancy
 * angle ψ (non-associated, ψ < φ).
 *
 * Solved by fem::nonsmooth_formulations::Plasticity with the numerics PlasticityProblem
 * (Drucker-Prager model, Voigt coordinates). See docs/devel/FEM_PLASTICITY_SCHEME.md.
 */
class DruckerPragerPlasticityNSL : public NonSmoothLaw {
 private:
  ACCEPT_SERIALIZATION(DruckerPragerPlasticityNSL);

  /** slope of the cone */
  double _eta{0.};

  /** dilatancy parameter (0: associated flow) */
  double _theta{0.};

  /** cohesion (in the cone coordinates) */
  double _cohesion{0.};

  /** default constructor */
  DruckerPragerPlasticityNSL() = default;

 public:
  /** constructor in the numerics parameters
   *
   *  \param eta slope of the cone (> 0)
   *  \param theta dilatancy parameter (0: associated flow)
   *  \param cohesion cohesion in the cone coordinates
   *  \param size size of the law (3: Voigt components of a 2D stress)
   */
  DruckerPragerPlasticityNSL(double eta, double theta, double cohesion, unsigned int size = 3)
      : NonSmoothLaw(size), _eta(eta), _theta(theta), _cohesion(cohesion) {};

  /** Mohr-Coulomb law of the plane (plane strain, σzz intermediate):
   *  R ≤ -p sin φ + c cos φ, dilatancy angle ψ
   *
   *  \param c cohesion
   *  \param phi friction angle (radians)
   *  \param psi dilatancy angle (radians, psi = phi: associated flow)
   *  \return the law, η = sin φ, θ = sin φ - sin ψ, cohesion √2 c cos φ
   */
  static std::shared_ptr<DruckerPragerPlasticityNSL> mohrCoulombPlaneStrain(double c,
                                                                            double phi,
                                                                            double psi);

  /** Destructor */
  ~DruckerPragerPlasticityNSL() noexcept = default;

  /** not implemented (throws) */
  bool isVerified(void) const override;

  /** \return the slope of the cone */
  inline double eta() const { return _eta; };

  /** \return the dilatancy parameter */
  inline double theta() const { return _theta; };

  /** \return the cohesion (cone coordinates) */
  inline double cohesion() const { return _cohesion; };

  /** print the data to the screen */
  void display() const override;

  // visitors hook
  virtual void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // DRUCKERPRAGERPLASTICITYNSL_H
