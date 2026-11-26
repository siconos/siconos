// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file  FEM-MoreauJeanOSI.hpp */

#ifndef FEMMoreauJeanOSI_H
#define FEMMoreauJeanOSI_H

#include "MoreauJeanOSI.hpp"

namespace siconos::mechanics::fem::integrators {

/**
   One Step time Integrator, Moreau-Jean algorithm dedicated to nsds based
   on DS and relations from the mechanics::fem component.

   For details regarding MoreauJeanOSI, see siconos::integrators::MoreauJeanOSI

*/
class MoreauJeanOSI : public siconos::integrators::MoreauJeanOSI {
 protected:
  ACCEPT_SERIALIZATION(MoreauJeanOSI);

 public:
  using siconos::integrators::MoreauJeanOSI::MoreauJeanOSI;

  /** This enum is used to get access to work vectors relared to DS
   *  It corresponds to:
   *  - a container saved in the graph of ds
   *  - a container associated to a specific ds
   */
  enum class wk_ds : std::size_t {
    residu_free,
    vfree,
    v_iter,
    residu_sigma_free,
    sigma_free,
    sigma_iter,
    q_sigma_free,
    buffer,
    size
  };

  /** destructor */
  virtual ~MoreauJeanOSI() noexcept = default;

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   *
   *  @param t time of initialization
   *  @param ds the dynamical system
   */
  virtual void initializeWorkVectorsForDS(
      double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   *
   *  \param inter the interaction
   *  \param interProp the properties on the graph
   *  \param DSG the dynamical systems graph
   */
  virtual void initializeWorkVectorsForInteraction(
      siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
      siconos::graphs::DynamicalSystemsGraph& DSG) override;

  /** initialize iteration matrix W MoreauJeanOSI matrix at time t
   *
   *  For a SolidLinearTIDS (unknowns \f$ x = [v; \sigma] \f$, compliance S, strain operator B)
   *  \f[ W = \begin{bmatrix} M & h\theta B^T \\ -h\theta B & S \end{bmatrix} \f]
   *  whose second row is the discretised constitutive law
   *  \f$ S\sigma_{k+1} - h\theta Bv_{k+1} = S\sigma_k + h(1-\theta)Bv_k - h\dot\epsilon_p \f$.
   *  W is not symmetric but its symmetric part diag(M, S) is positive definite, so the
   *  Delassus operator of the stress interactions is positive (with the symmetric
   *  indefinite form [M, h theta B^T; h theta B, -S] it is negative and the plasticity
   *  problem has no physical solution). See docs/devel/FEM_PLASTICITY_SCHEME.md.
   *
   *  @param time
   *  @param ds a pointer to DynamicalSystem
   */
  /** \return true if the LU factorization of the iteration matrix is needed (computeFreeState,
   *  integrate and computeIteration solve with it); false for the global integrator, where
   *  the one-step problem holds W */
  virtual bool needsIterationMatrixLU() const { return true; }

  virtual void initializeIterationMatrix(
      double time, std::shared_ptr<siconos::modeling::SecondOrderDS> ds) const override;

  /** @return the inverse of the iteration matrix
   *  @param ds the concerned dynamical system
   *  @param LUW the LU-factorisation of W
   */
  std::shared_ptr<siconos::algebra::SiconosMatrix> iterationMatrixInverse(
      std::shared_ptr<siconos::modeling::SecondOrderDS> ds,
      const siconos::algebra::SiconosDenseLUMatrix& LUW);

  /** compute the initial state of the Newton loop */
  virtual void computeInitialNewtonState() override;

  /**
      return the maximum of all norms for the "MoreauJeanOSI-discretized"
      residus of DS

      @return a double
    */
  double computeResidu() override;

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  virtual void computeFreeState() override;

  /** integrate the system, between tinit and tend (->iout=true), with possible
   *  stop at tout (->iout=false)
   *
   *  @param tinit the initial time
   *  @param tend the end time
   *  @param tout the real end time
   *  @param notUsed useless flag (for MoreauJeanOSI, used in LsodarOSI)
   */
  void integrate(double& tinit, double& tend, double& tout, int& notUsed) override;

  /** update the state of the dynamical systems
   *
   *  @param level the level of interest for the dynamics: not used at the time
   */
  virtual void updateState(const unsigned int level) override;

  /** compute the nonsmooth inputs from the multipliers (see MoreauJeanOSI::updateInput),
   *  then convert, for each SolidLinearTIDS, the input of the stress interactions
   *  (\f$ C^T\lambda \f$, accumulated in the plastic rate vector) into the plastic
   *  deformation rate \f$ \dot\epsilon_p = -C^T\lambda / h \f$.
   *
   *  The stress rows of the iteration matrix (see initializeIterationMatrix) read
   *  \f$ S\sigma_{k+1} - h\theta B v_{k+1} = S\sigma_k + h(1-\theta)Bv_k - h\dot\epsilon_p \f$
   *  and the multiplier enters them as \f$ Hr = C^T\lambda \f$. With an interaction
   *  \f$ y = C\sigma + c \geq 0 \f$ (yield function \f$ f = -y \f$), this is the associated
   *  flow rule \f$ \dot\epsilon_p = \lambda/h\, \partial f/\partial\sigma \f$.
   *  The simulation resets the plastic rate before each call (see
   *  SolidLinearTIDS::resetNonSmoothPart).
   *
   *  @param time current time
   *  @param level level of the input
   */
  void updateInput(double time, unsigned int level) override;
  using siconos::integrators::MoreauJeanOSI::updateInput;

  /** compute the current iteration
   *
   */
  void computeIteration() override;

  /** Apply the rule to one Interaction to know if it should be included in the
   *  IndexSet of level i
   *
   *  @param inter the Interaction to test
   *  @param i level of the IndexSet
   *  @return Boolean
   */
  virtual bool addInteractionInIndexSet(
      std::shared_ptr<siconos::modeling::Interaction> inter,
      siconos::graphs::InteractionsGraph::size_type i) override;
};

}  // namespace siconos::mechanics::fem::integrators
#endif  // MoreauJeanOSI_H
