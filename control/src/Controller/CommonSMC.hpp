/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
  \brief General interface to define a sliding mode actuator
  */

#ifndef CommonSMC_H
#define CommonSMC_H

#include "Actuator.hpp"
#include "SiconosAlgebraTypeDef.hpp"

#include <relay_cst.h>

class CommonSMC : public Actuator
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(CommonSMC);

protected:
  /** default constructor */
  CommonSMC() {};

  /** index for saving data */
  unsigned int _indx;

  /** name of the plugin to add a term to the sliding variable; useful when doing trajectory tracking */
  std::string _plugineName;

  /** name of the plugin to compute \f$y = h(x, ...)\f$ for the nonlinear case*/
  std::string _pluginhName;

  /** name of the plugin to compute \f$\nabla_x h\f$ for the nonlinear case*/
  std::string _pluginJachxName;

  /** name of the plugin to compute \f$\nabla_\lambda h\f$ for the nonlinear case*/
  std::string _pluginJachlambdaName;

  /** name of the plugin to compute \f$\nabla_\lambda g\f$ for the nonlinear case*/
  std::string _pluginJacglambdaName;

  /** the vector defining the linear contribution of the state to the sliding variable  (\f$ \sigma = Cx \f$) */
  SP::SimpleMatrix _Csurface;

  /** matrix describing the influence of \f$lambda\f$ on \f$\sigma\f$ */
  SP::SimpleMatrix _D;

  /** scalar multiplying Sign; \f$ u^s = - \alpha Sign \f$ */
  double _alpha;

  /** the Relation for the Controller */
  SP::FirstOrderR _relationSMC;

  /** Interaction for the control */
  SP::Interaction _interactionSMC;

  /** easy access to lambda */
  SP::SiconosVector _lambda;

  /** Relay solver type */
  int _numericsSolverId;

  /** Numerical precision expected for the Relay solver */
  double _precision;

  /** the Model for the controller */
  SP::Model _SMC;

  /** the DynamicalSystem for the controller */
  SP::FirstOrderNonLinearDS _DS_SMC; // XXX replace this by FirstOrderDS

  /** the TimeDiscretisation for the controller */
  SP::TimeDiscretisation _td;

  /** Simulation for the controller */
  SP::TimeStepping _simulationSMC;

  /** Integrator for the controller */
  SP::OneStepIntegrator _integratorSMC;

  /** Theta for the controller */
  double _thetaSMC;

  /** OneStepNsProblem for the controller */
  SP::LinearOSNS _OSNSPB_SMC;

  /** SP::EventsManager of the SMC Simulation */
  SP::EventsManager _eventsManager;

  /** SP::NonSmoothLaw for computing the control law */
  SP::NonSmoothLaw _nsLawSMC;

  /** inverse of CB */
  SP::SimpleMatrix _invCB;

  /** Store \f$u^{eq}\f$ */
  SP::SiconosVector _ueq;

  /** Store \f$u^s\f$ */
  SP::SiconosVector _us;

  /** Do not use the state-continuous equivaluent control \f$u^{eq}\f$ */
  bool _noUeq;

  /** If true perform the computation of the residus in the Newton loop if needed */
  bool _computeResidus;

  /** Compute the equivalent part of the control \f$u^{eq}\f$.
   * The method used here is to discretize the continuous-time
   * formula using a theta method
   */
  void computeUeq();

public:

  /** Constructor with a TimeDiscretisation and a Model.
   * \param type the type of the SMC Actuator
   * \param sensor the ControlSensor feeding the Actuator
   */
  CommonSMC(unsigned int type, SP::ControlSensor sensor): Actuator(type, sensor),
    _indx(0), _alpha(1.0), _numericsSolverId(SICONOS_RELAY_AVI_CAOFERRIS), _precision(1e-8),
    _thetaSMC(0.5), _noUeq(false), _computeResidus(true) {}

  /** Constructor with a TimeDiscretisation, a Model and two matrices
   * \param type the type of the SMC Actuator
   * \param sensor the ControlSensor feeding the Actuator
   * \param B the B matrix
   * \param D the saturation matrix
   */
  CommonSMC(unsigned int type, SP::ControlSensor sensor, SP::SimpleMatrix B, SP::SimpleMatrix D = std11::shared_ptr<SimpleMatrix>()):
    Actuator(type, sensor), _indx(0), _D(D), _alpha(1.0), _numericsSolverId(SICONOS_RELAY_AVI_CAOFERRIS),
    _precision(1e-8), _thetaSMC(0.5), _noUeq(false), _computeResidus(true)
  {
    _B = B;
  }


  /** Compute the new control law at each event
  */
  virtual void actuate() = 0;

  /** Initialization
   * \param m the Model
   */
  virtual void initialize(const Model& m);


  void sete(const std::string& plugin);
  void seth(const std::string& plugin);
  void setJachx(const std::string& plugin);
  void setJachlambda(const std::string& plugin);
  void setg(const std::string& plugin);
  void setJacgx(const std::string& plugin);
  void setJacglambda(const std::string& plugin);

  /** Set Csurface
   * \param Csurface a SP::SimpleMatrix containing the new value for _Csurface
   */
  void setCsurface(SP::SimpleMatrix Csurface);

  /** Set _D to pointer newPtr
   * \param newSat a SP::SimpleMatrix containing the new value for _D
   */
  void setSaturationMatrix(SP::SimpleMatrix newSat);

  /** Set _alpha
   * \param alpha the new value for _alpha
   */
  inline void setAlpha(double alpha) { _alpha = alpha; };

  /** get _lambda
   * \return a pointer to _lambda
   */
  inline SP::SiconosVector lambda() const
  {
    return _lambda;
  };

  /** Set the solver
   * \param numericsSolverId the solver for the relay
   */
  inline void setSolver(const int numericsSolverId)
  {
    _numericsSolverId = numericsSolverId;
  };

  /** Set the precision
   * \param newPrecision a double
   */
  inline void setPrecision(double newPrecision)
  {
    _precision = newPrecision;
  };

  /** Get the OneStepNSProblem problem associated with the controller. This is useful to
   * gain access to the data given to the solver in Numerics
   * \return a reference to the LinearOSNS problem
   */
  inline const LinearOSNS& relay()
  {
    return * _OSNSPB_SMC;
  };

  /** get \f$u^{eq}\f$
   * \return a reference to _ueq
   */
  inline SiconosVector& ueq()
  {
    return *_ueq;
  };

  /** get \f$u^{s}\f$
   * \return a reference to _us
   */

  inline SiconosVector& us()
  {
    return *_us;
  };

  /** Set _theta, used in some discretization method for \f$u^{eq}\f$
   * \param newTheta the new value for _thetaSMC
   */

  inline void setTheta(double newTheta)
  {
    _thetaSMC = newTheta;
  };

  /** Disable (or enable) the use of the state-continuous control \f$u^{eq}\f$
   * \param b disable the use of Ueq if true
   */
  inline void noUeq(bool b)
  {
    _noUeq = b;
  };

  /** Disable (or enable) the computation of the residus on the Newton loop.
   * This has an incidence only if the Relation is nonlinear
   * \param b disable the computation of the residus
   */
  inline void setComputeResidus(bool b)
  {
    _computeResidus = b;
  };

  /** This is derived in child classes if they need to copy the TimeDiscretisation
   * associated with this Sensor
  *  \param td the TimeDiscretisation for this Sensor
  */
  virtual void setTimeDiscretisation(const TimeDiscretisation& td);

  /** Set the DynamicalSystem used to compute the control law.
   * This is useful when we have a Nonlinear problem and we need to compute
   * the control law with an approximate model, or when the dynamics are
   * quite different.
   * \param ds the DynamicalSystem to be used in the Controller
   */
  void setDS(SP::FirstOrderNonLinearDS ds) // XXX replace this by FirstOrderDS
  {
    _DS_SMC = ds;
  };

  /** get the Model used in the SMC
   * \return the Model used in the SMC
   */
  virtual SP::Model getInternalModel() const { return _SMC; };

  /** get the Integrator used in the SMC
   * \return the Integrator used in the SMC
   */
  OneStepIntegrator& getInternalOSI() const { return *_integratorSMC; };

};
#endif
