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

/*! \file CommonSMC.hpp
  General interface to define a sliding mode actuator
*/

#ifndef CommonSMC_H
#define CommonSMC_H

#include "Actuator.hpp"
#include "FunctionTypes.hpp"
#include "relay_cst.h"  // contains only enum. Ok.

namespace siconos::modeling {

class FirstOrderR;
class Interaction;
class FirstOrderNonLinearDS;
class NonSmoothLaw;
}  // namespace siconos::modeling

namespace siconos::integrators {
class OneStepIntegrator;
}

namespace siconos::nonsmooth_formulations {
class LinearOSNS;
}  // namespace siconos::nonsmooth_formulations

namespace siconos::simulation {
class TimeStepping;
class EventsManager;

}  // namespace siconos::simulation
namespace siconos::control {

class CommonSMC : public Actuator {
 private:
  ACCEPT_SERIALIZATION(CommonSMC);

 protected:
  /** index for saving data */
  unsigned int _indx{0};

  /** function wrapper used to compute e(t) in FirstOrderLinearR */
  siconos::modeling::func_prototypes::FunctionS_V computee_{nullptr};

  /** function wrapper used to compute \f$ h(x,t,\lambda) in relation \f$ */
  siconos::modeling::func_prototypes::FunctionBVSV_V computeh_{nullptr};

  /** function wrapper used to compute  \f$ \nabla_x h(x,t,\lambda)  in relation\f$ */
  siconos::modeling::func_prototypes::FunctionBVSV_M computejacobianhOver_state_{nullptr};

  /** function wrapper used to compute  \f$ \nabla_{\lambda} h(x,t,\lambda)  in relation\f$ */
  siconos::modeling::func_prototypes::FunctionBVSV_M computejacobianhOver_lambda_{nullptr};

  /** function wrapper used to compute  \f$ \nabla_{\lambda} g(x,\lambda)  in relation\f$ */
  siconos::modeling::func_prototypes::FunctionBVSV_M computejacobiangOver_lambda_{nullptr};

  // Note: g and  \f$ \nabla_x g(x,\lambda) \f$ */ are attributes of the Actuator base class.

  /** the vector defining the linear contribution of the state to the sliding variable  ( \f$
   * \sigma = Cx \f$ ) */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _Csurface{nullptr};

  /** matrix describing the influence of \f$ lambda \f$  on \f$ \sigma \f$ */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _D{nullptr};

  /** scalar multiplying Sign; \f$ u^s = - \alpha Sign \f$ */
  double _alpha{1.};

  /** the Relation for the Controller */
  std::shared_ptr<siconos::modeling::FirstOrderR> _relationSMC{nullptr};

  /** Interaction for the control */
  std::shared_ptr<siconos::modeling::Interaction> _interactionSMC{nullptr};

  /** easy access to lambda */
  std::shared_ptr<siconos::algebra::SiconosVector> _lambda{nullptr};

  /** Relay solver type */
  int _numericsSolverId{SICONOS_RELAY_AVI_CAOFERRIS};

  /** Numerical precision expected for the Relay solver */
  double _precision{1e-8};

  /** the nsds for the controller */
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _nsdsSMC{nullptr};

  /** the DynamicalSystem for the controller */  // XXX replace this by FirstOrderDS
  std::shared_ptr<siconos::modeling::FirstOrderNonLinearDS> _DS_SMC{nullptr};

  /** Internal buffer for b vector of DS_SMC_ */
  siconos::algebra::SiconosVector bSMC_;

  /** the TimeDiscretisation for the controller */
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _td{nullptr};

  /** Simulation for the controller */
  std::shared_ptr<siconos::simulation::TimeStepping> _simulationSMC{nullptr};

  /** Integrator for the controller */
  std::shared_ptr<siconos::integrators::OneStepIntegrator> _integratorSMC{nullptr};

  /** Theta for the controller */
  double _thetaSMC{0.5};

  /** OneStepNsProblem for the controller */
  std::shared_ptr<siconos::nonsmooth_formulations::LinearOSNS> _OSNSPB_SMC{nullptr};

  /** std::shared_ptr<siconos::simulation::EventsManager> of the SMC Simulation */
  std::shared_ptr<siconos::simulation::EventsManager> _eventsManager{nullptr};

  /** std::shared_ptr<siconos::modeling::NonSmoothLaw> for computing the control law */
  std::shared_ptr<siconos::modeling::NonSmoothLaw> _nsLawSMC{nullptr};

  /** inverse of CB */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _invCB{nullptr};

  /** Store  \f$ u^{eq} \f$  */
  std::shared_ptr<siconos::algebra::SiconosVector> _ueq{nullptr};

  /** Store  \f$ u^s \f$  */
  std::shared_ptr<siconos::algebra::SiconosVector> _us{nullptr};

  /** Do not use the state-continuous equivaluent control  \f$  u^{eq} \f$ */
  bool _noUeq{false};

  /** If true perform the computation of the residus in the Newton loop if needed */
  bool _computeResidus{true};

  /** Compute the equivalent part of the control  \f$  u^{eq} \f$.
   *  The method used here is to discretize the continuous-time
   *  formula using a theta method
   */
  void computeUeq();

 public:
  /** General constructor
   *
   *  \param type the type of the SMC Actuator
   *  \param sensor the ControlSensor feeding the Actuator
   */
  CommonSMC(ActuatorType type, std::shared_ptr<ControlSensor> sensor)
      : Actuator(type, sensor) {}

  /** Constructor for dynamics affine in control
   *
   *  \param type the type of the SMC Actuator
   *  \param sensor the ControlSensor feeding the Actuator
   *  \param B the matrix multiplying the control input
   *  \param D the saturation matrix (optional)
   */
  CommonSMC(ActuatorType type, std::shared_ptr<ControlSensor> sensor,
            std::shared_ptr<siconos::algebra::SiconosMatrix> B,
            std::shared_ptr<siconos::algebra::SiconosMatrix> D = nullptr)
      : Actuator(type, sensor, B), _D(D) {}

  /** Compute the new control law at each event
   */
  virtual void actuate() = 0;

  /** Initialization
   *
   *  \param nsds current nonsmooth dynamical system
   *  \param s current simulation setup
   */
  virtual void initialize(const siconos::modeling::NonSmoothDynamicalSystem& nsds,
                          const siconos::simulation::Simulation& s);

  void sete(const siconos::modeling::func_prototypes::FunctionS_V&
                fct);  // Meaningful only for FirstOrderLinearR

  /** set a user-defined function to compute \f$ h(x,t,\lambda) \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputehFunction(const siconos::modeling::func_prototypes::FunctionBVSV_V& fct);

  /** set a user-defined function to compute \f$ \nabla_x h(x, t, \lambda) \f$ \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianhOver_stateFunction(
      const siconos::modeling::func_prototypes::FunctionBVSV_M& fct);

  /** set a user-defined function to compute \f$ \nabla_{\lambda} h(x, t, \lambda) \f$ \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobianhOver_lambdaFunction(
      const siconos::modeling::func_prototypes::FunctionBVSV_M& fct);

  /** set a user-defined function to compute \f$ \nabla_{\lambda} g(x, \lambda) \f$ \f$
   *
   *  \param fct the user-defined function (std::function, lambda ...)
   */
  void setComputeJacobiangOver_lambdaFunction(
      const siconos::modeling::func_prototypes::FunctionBVSV_M& fct);

  /** Set Csurface
   *
   *  \param Csurface a std::shared_ptr<siconos::algebra::SiconosMatrix> containing the new
   * value for _Csurface
   */
  void setCsurface(std::shared_ptr<siconos::algebra::SiconosMatrix> Csurface);

  /** Set _D to pointer newPtr
   *
   *  \param newSat a std::shared_ptr<siconos::algebra::SiconosMatrix> containing the new value
   * for _D
   */
  void setSaturationMatrix(std::shared_ptr<siconos::algebra::SiconosMatrix> newSat);

  /** Set _alpha
   *
   *  \param alpha the new value for _alpha
   */
  inline void setAlpha(double alpha) { _alpha = alpha; };

  /** get _lambda
   *
   *  \return a pointer to _lambda
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> lambda() const { return _lambda; };

  /** Set the solver
   *
   *  \param numericsSolverId the solver for the relay
   */
  inline void setSolver(const int numericsSolverId) { _numericsSolverId = numericsSolverId; };

  /** Set the precision
   *
   *  \param newPrecision a double
   */
  inline void setPrecision(double newPrecision) { _precision = newPrecision; };

  /** Get the OneStepNSProblem problem associated with the controller. This is useful to
   *  gain access to the data given to the solver in Numerics
   *
   *  \return a reference to the LinearOSNS problem
   */
  inline const siconos::nonsmooth_formulations::LinearOSNS& relay() { return *_OSNSPB_SMC; };

  /** get \f$ u^{eq} \f$
   *
   *  \return a reference to _ueq
   */
  inline siconos::algebra::SiconosVector& ueq() { return *_ueq; };

  /** get  \f$ u^{s} \f$
   *
   *  \return a reference to _us
   */

  inline siconos::algebra::SiconosVector& us() { return *_us; };

  /** Set _theta, used in some discretization method for \f$ u^{eq} \f$
   *
   *  \param newTheta the new value for _thetaSMC
   */

  inline void setTheta(double newTheta) { _thetaSMC = newTheta; };

  /** Disable (or enable) the use of the state-continuous control \f$ u^{eq} \f$
   *
   *  \param b disable the use of Ueq if true
   */
  inline void noUeq(bool b) { _noUeq = b; };

  /** Disable (or enable) the computation of the residus on the Newton loop.
   *  This has an incidence only if the Relation is nonlinear
   *
   *  \param b disable the computation of the residus
   */
  inline void setComputeResidus(bool b) { _computeResidus = b; };

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   *  associated with this Sensor
   *
   *  \param td the TimeDiscretisation for this Sensor
   */
  virtual void setTimeDiscretisation(const siconos::simulation::TimeDiscretisation& td);

  /** Set the DynamicalSystem used to compute the control law.
   *  This is useful when we have a Nonlinear problem and we need to compute
   *  the control law with an approximate model, or when the dynamics are
   *  quite different.
   *
   *  \param ds the DynamicalSystem to be used in the Controller
   */
  void setDS(std::shared_ptr<siconos::modeling::FirstOrderNonLinearDS>
                 ds)  // XXX replace this by FirstOrderDS
  {
    _DS_SMC = ds;
  };

  /** get the NSDS used in the SMC
   *
   *  \return the NSDS used in the SMC
   */
  virtual std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> getInternalNSDS()
      const {
    return _nsdsSMC;
  };

  /** get the Integrator used in the SMC
   *
   *  \return the Integrator used in the SMC
   */
  siconos::integrators::OneStepIntegrator& getInternalOSI() const { return *_integratorSMC; };
};

}  // namespace siconos::control
#endif
