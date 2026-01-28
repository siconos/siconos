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
/*! \file
  NewMark Alpha Scheme Time-Integrator for Dynamical Systems
*/
#ifndef NEWMARKALPHAOSI_H
#define NEWMARKALPHAOSI_H

#include "OneStepIntegrator.hpp"

namespace siconos::modeling {

class LagrangianDS;

}
namespace siconos::integrators {

/////// Free Functions ///////
namespace newmark_alpha {

// Names of the Newmark integrator parameters
// Used for newmark array indices;
namespace names {
constexpr int Beta = 0;
constexpr int Gamma = 1;
constexpr int Alpha_m = 2;
constexpr int Alpha_f = 3;
}  // namespace names

/** compute NewmarkAlpha W iteration matrix for Lagrangian systems
 *
 *  \param[in] time current time
 *  \param[in] h current time-step
 *  \param[in] params integrator parameters ([beta, gamma, alpha_m, alpha_f])
 *  \param[in] ds a Lagrangian dynamical system
 *  \param[in,out] W the result in W
 *  \param[in,out] LUW the LU factorisation of W (updated)
 */
void computeIterationMatrix_Lagrangian(
    double time, double h, const std::array<double, 4>& params,
    siconos::modeling::LagrangianDS& ds, Eigen::Ref<siconos::algebra::SiconosMatrix> W,
    std::shared_ptr<siconos::algebra::SiconosDenseLUMatrix>& LUW);

}  // namespace newmark_alpha

/**  NewMarkAlpha Scheme Time-Integrator for Dynamical Systems
 *
 *
 * NewMarkAlphaOSI is used to solve constrained dynamical systems represented by
 * index-3 DAE
 *
 * NewMarkAlphaOSI is instantiated with values of beta, gamma, alpha_m, alpha_f
 * and the list of concerned dynamical systems. Each DynamicalSystem is
 * associated to a SiconosMatrix named "W"
 *
 * W matrices are initialized and computed in initializeIterationMatrix and
 * computeIterationMatrix.
 */
class NewMarkAlphaOSI : public OneStepIntegrator {
 protected:
  ACCEPT_SERIALIZATION(NewMarkAlphaOSI);
  /** Parameters of the numerical scheme:  beta, gamma, alpha_m, alpha_f */
  std::array<double, 4> newmark_{};  // all 0. by default.

  /** Order of the polynomial for dense output*/
  unsigned int denseOutputPolynomialOrder_{5};

  /** Indicator whether or not constraints at the velocity level are handled
   *  true: constraints at the velocity level are handled
   *  false: constraints at the position are handled
   */
  bool handleVelocityConstraints_{false};

 public:
  enum NewMarkAlphaOSI_ds_workVector_id {
    RESIDU_FREE,
    FREE,
    ACCE_LIKE,
    ACCE_MEMORY,
    WORK_LENGTH
  };

  enum NewMarkAlphaOSI_interaction_workVector_id { OSNSP_RHS, WORK_INTERACTION_LENGTH };

  enum NewMarkAlphaOSI_workBlockVector_id { xfree, BLOCK_WORK_LENGTH };

  enum NewMarkAlphaOSI_interaction_workMat_id { DENSE_OUTPUT_COEFFICIENTS, MAT_WORK_LENGTH };

  /** constructor with only parameters beta, gamma, alpha_m, alpha_f
   * \param beta double
   * \param gamma double
   * \param alpha_m double
   * \param alpha_f double
   * \param flag true of working at velocity level
   */
  NewMarkAlphaOSI(double beta, double gamma, double alpha_m, double alpha_f, bool flag);

  /** constructor with only the parameter rho_infty
   * \param rho_infty double
   * \param flag true of working at velocity level
   */
  NewMarkAlphaOSI(double rho_infty, bool flag);

  /** destructor
   */
  virtual ~NewMarkAlphaOSI() noexcept = default;

  // --- GETTERS/SETTERS ---

  /** set value to the parameter beta
   * \param beta value of beta
   */
  inline void setBeta(double beta) { newmark_[newmark_alpha::names::Beta] = beta; };

  /** set value to the parameter gamma
   * \param value_gamma double : value of gamma
   */
  inline void setGamma(double value_gamma) {
    newmark_[newmark_alpha::names::Gamma] = value_gamma;
  };

  /** set value to the parameter alpha_m
   * \param value_alpha_m double : value of alpha_m
   */

  inline void setAlpha_m(double value_alpha_m) {
    newmark_[newmark_alpha::names::Alpha_m] = value_alpha_m;
  };

  /** set value to the parameter alpha_f
   * \param value_alpha_f double : value of alpha_f
   */

  inline void setAlpha_f(double value_alpha_f) {
    newmark_[newmark_alpha::names::Alpha_f] = value_alpha_f;
  };

  /** set values to the parameters beta, gamma, alpha_f, alpha_m from the value
   * of rho_infty \param rho_infty double : value of rho_infty
   */

  inline void setParametersFromRho_infty(double rho_infty) {
    newmark_[newmark_alpha::names::Alpha_m] = (2 * rho_infty - 1) / (rho_infty + 1);
    newmark_[newmark_alpha::names::Alpha_f] = rho_infty / (rho_infty + 1);
    newmark_[newmark_alpha::names::Gamma] = 0.5 + newmark_[newmark_alpha::names::Alpha_f] -
                                            newmark_[newmark_alpha::names::Alpha_m];
    newmark_[newmark_alpha::names::Beta] =
        0.25 * std::pow((newmark_[newmark_alpha::names::Gamma] + 0.5), 2);
  };

  /** get value of beta
   * \return double
   */
  inline double getBeta() { return newmark_[newmark_alpha::names::Beta]; };

  /** get value of gamma
   * \return double
   */
  inline double getGamma() { return newmark_[newmark_alpha::names::Gamma]; };

  /** get value of alpha_m
   * \return double
   */
  inline double getAlpha_m() { return newmark_[newmark_alpha::names::Alpha_m]; };

  /** get value of alpha_f
   * \return double
   */

  inline double getAlpha_f() { return newmark_[newmark_alpha::names::Alpha_f]; };

  /** \return the order of the polynomial for dense output
   */
  inline unsigned int denseOutputPolynomialOrder() const {
    return denseOutputPolynomialOrder_;
  }

  /** set to true if constraints at the velocity level are handled, else false (default)
   * \param flag bool
   */
  inline void setHandleVelocityConstraints(bool flag) { handleVelocityConstraints_ = flag; }

  /** get the flag handleVelocityConstraints_
   * \return bool
   */
  inline bool handleVelocityConstraints() const { return handleVelocityConstraints_; }

  /** initialize W matrix
   *  \param time current time value
   *  \param ds a pointer to DynamicalSystem
   */
  void initializeIterationMatrix(double time,
                                 std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** compute the residual of dynamical equation
   *\return double: maximum residu over all DSs
   */
  double computeResidu() override;

  /** compute the free state of the discretized dynamical system */
  void computeFreeState() override;

  /** integrates the Interaction linked to this integrator, without taking
   * non-smooth effects into account \param vertex_inter of the interaction
   * graph \param osnsp pointer to OneStepNSProblem
   */
  void computeFreeOutput(siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
                         siconos::nonsmooth_formulations::OneStepNSProblem* osnsp) override;

  /** return the workVector corresponding to the right hand side of the OneStepNonsmooth
   * problem
   */
  siconos::algebra::SiconosVector& osnsp_rhs(
      siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
      siconos::graphs::InteractionsGraph& indexSet) override {
    return *(*indexSet.properties(vertex_inter).workVectors)[NewMarkAlphaOSI::OSNSP_RHS];
  };

  /** initialize */
  //  void initialize(Model& m);

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   * \param t time of initialization
   * \param ds the dynamical system
   */
  void initializeWorkVectorsForDS(
      double t, std::shared_ptr<siconos::modeling::DynamicalSystem> ds) override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   * \param inter the interaction
   * \param interProp the properties on the graph
   * \param DSG the dynamical systems graph
   */
  void initializeWorkVectorsForInteraction(
      siconos::modeling::Interaction& inter, siconos::graphs::InteractionProperties& interProp,
      siconos::graphs::DynamicalSystemsGraph& DSG) override;

  /** get the number of index sets required for the simulation
   * \return unsigned int
   */
  unsigned int numberOfIndexSets() const override { return 3; };

  /** prepare for Newton Iteration
   * \param time
   */
  void prepareNewtonIteration(double time) override;

  /** predict first values for the Newton iteration */
  void prediction();

  /** correct state of all levels of Dynamical Systems after each Newton
   * iteration
   */
  void correction();

  /** integrate the system, between tinit and tend with possible stop at tout
   *  \param tinit double: tinit, initial time
   *  \param tend double: tend, end time
   *  \param tout double: tout, real end time
   *  \param flag useless for NewMarkAlphaOSI
   */
  void integrate(double& tinit, double& tend, double& tout, int& flag) override;

  /** updates the state of the Dynamical Systems
   *  \param level the level of interest for the dynamics: not used at the time
   */
  void updateState(const unsigned int level) override;

  /** compute the current iteration
   *
   */
  void computeIteration() override { assert(0); };
  /** Compute coefficients of the polynomial of the dense output for a given DS
   *  \param ds std::shared_ptr<siconos::modeling::DynamicalSystem>, ds concerned
   */
  void computeCoefsDenseOutput(std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /** prepare for Event localization*/
  void prepareEventLocalization();

  /** Generate dense output for all Dynamical Systems belonging to OSI
   * \param time at which we want to generate the dense output
   */
  void DenseOutputallDSs(double time);

  /** Displays the data of the NewMarkAlpha's integrator
   */
  void display() const override;
};

}  // namespace siconos::integrators
#endif  // NEWMARKALPHAOSI_H
