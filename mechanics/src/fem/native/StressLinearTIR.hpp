// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file StressLinearTIR.hpp

 */
#ifndef STRESSLINEARRELATION_H
#define STRESSLINEARRELATION_H

#include "LagrangianLinearTIR.hpp"

namespace siconos::mechanics::fem {
/**
    Stress Linear Relation.

    Stress Relation with:

    \f$ y= [0 C] [v ; \sigma] + e\f$

    \f$ p = [0 C]^t \lambda = [0 C^t \lambda] \f$

    With \f$ y \geq 0 \f$ (yield function \f$ f = -y \f$), \f$ C^T\lambda \f$ is the
    plastic increment that enters the stress rows of the iteration matrix
    (see fem::integrators::MoreauJeanOSI::initializeIterationMatrix).
    computeInput accumulates it in SolidLinearTIDS::plasticRate() (several interactions
    may share a solid; the solid resets it before each input update), and
    fem::integrators::MoreauJeanOSI::updateInput converts it into the plastic deformation
    rate \f$ \dot\epsilon_p = -C^T\lambda / h \f$.
    See docs/devel/FEM_PLASTICITY_SCHEME.md.

*/
class StressLinearTIR : public siconos::modeling::LagrangianLinearTIR {
 protected:
  ACCEPT_SERIALIZATION(StressLinearTIR);

 public:
  /** Names and positions of ds variables that might be used (read only) in compute input and
   * output.
   * Last param, size, gives the size of the BlockVector to be allocated.
   */
  enum class ds_var : std::size_t {
    sigma,
    sigma1,
    sigma2,
    epsilon_p,
    epsilon_p1,
    epsilon_p2,
    size
  };

  using LagrangianLinearTIR::LagrangianLinearTIR;

  StressLinearTIR(Eigen::Ref<siconos::algebra::SiconosMatrix> C);

  /** destructor
   */
  virtual ~StressLinearTIR() noexcept = default;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an siconos::modeling::Interaction using this relation
   */
  void checkSize(const siconos::modeling::Interaction& inter) const override;

  /** default function to compute y
   *
   *  \param time not used
   *  \param inter the siconos::modeling::Interaction we want to update
   *  \param derivativeNumber the derivative of y we want to compute
   */
  void computeOutput(double time, siconos::modeling::Interaction& inter,
                     siconos::algebra::blocks::size_type derivativeNumber = 0) override;

  /** default function to compute r
   *
   *  \param time not used
   *  \param inter the siconos::modeling::Interaction we want to update
   *  \param level the derivative of lambda we want to compute
   */
  void computeInput(double time, siconos::modeling::Interaction& inter,
                    siconos::algebra::blocks::size_type level = 0) override;

  /** compute all the H Jacobian
   *
   *  \param time not used
   *  \param inter the siconos::modeling::Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJach(double time, siconos::modeling::Interaction& inter) override {}

  /** compute all the G Jacobian
   *
   *  \param time not used
   *  \param inter the siconos::modeling::Interaction we want to update
   *  \param interProp interaction properties
   */
  void computeJacg(double time, siconos::modeling::Interaction& inter) override {}

  /** print the data to the screen
   */
  void display() const override;

  /** @brief allocation of memory space for relations in the graph
   *      Warning: internal use only (called from Topology)
   *  @param dslink a container of vectors (pointers), from the parent interaction
   *  @param ds1 first ds concerned by the relation
   *  @param ds2 second ds concerned by the relation
   */
  virtual void allocate_read_dynamical_systems_var_vectors(
      std::vector<std::shared_ptr<siconos::algebra::BlockVector>>& ds_vars,
      modeling::DynamicalSystem& ds1, modeling::DynamicalSystem& ds2) const override;
};

}  // namespace siconos::mechanics::fem
#endif
