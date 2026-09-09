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

/*! \file ControlSimulation.hpp
  \brief Abstract class - General interface for all Control Dynamical Systems.
*/

#ifndef CONTROLDYNAMICALSYSTEM_H
#define CONTROLDYNAMICALSYSTEM_H

#include <memory>
#include <string>

#include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"

namespace siconos::graphs {
struct DynamicalSystemsGraph;
struct InteractionsGraph;
}  // namespace siconos::graphs
namespace siconos::modeling {
class NonSmoothDynamicalSystem;
class DynamicalSystem;
}  // namespace siconos::modeling

namespace siconos::integrators {
class OneStepIntegrator;
}

namespace siconos::simulation {
class TimeDiscretisation;
class Simulation;
}  // namespace siconos::simulation

namespace siconos::control {

class ControlManager;
class Sensor;
class Actuator;
class Observer;

class ControlSimulation {
 private:
  ACCEPT_SERIALIZATION(ControlSimulation);

 protected:
  /** Constructor with the minimal set of data
   * \param t0 the starting time \f$ t_0 \f$
   * \param T the end time T
   * \param h the simulation time step
   * */
  ControlSimulation(double t0, double T, double h);

  /** Starting time */
  double _t0{0.};

  /** End time */
  double _T{0.};

  /** Simulation step */
  double _h{0.};

  /** Theta for MoreauJeanOSI */
  double _theta{0.5};

  /** Time spent computing */
  double _elapsedTime{0.};

  /** rough estimation of the number of points to save */
  size_t _N{0};

  /** Dimension of the state space */
  size_t _nDim{0};

  /** Save only the data in the main Simulation*/
  bool _saveOnlyMainSimulation{false};

  /** If true, do not show progress of the simulation */
  bool _silent{false};

  /** Matrix for saving result */
  std::shared_ptr<siconos::algebra::SiconosMatrix> _dataM{nullptr};

  /** Legend for the columns in the matrix _dataM*/
  std::string _dataLegend{"none"};

  /** NonSmoothDynamicalSystem */
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _nsds{nullptr};

  /** TimeDiscretisation for the simulation*/
  std::shared_ptr<siconos::simulation::TimeDiscretisation> _processTD{nullptr};

  /** The Simulation */
  std::shared_ptr<siconos::simulation::Simulation> _processSimulation{nullptr};

  /** The integrator */
  std::shared_ptr<siconos::integrators::OneStepIntegrator> _processIntegrator{nullptr};

  /** the ControlManager */
  std::shared_ptr<ControlManager> _CM{nullptr};

  /** siconos::graphs::DynamicalSystemsGraph (for convenience)*/
  std::shared_ptr<siconos::graphs::DynamicalSystemsGraph> _DSG0{nullptr};

  /** InteractionsGraph (for convenience)*/
  std::shared_ptr<siconos::graphs::InteractionsGraph> _IG0{nullptr};

  // Rule of five
  ControlSimulation() = delete;
  ControlSimulation(const ControlSimulation&) = delete;
  ControlSimulation(ControlSimulation&&) = delete;
  ControlSimulation& operator=(const ControlSimulation&) = delete;
  ControlSimulation& operator=(ControlSimulation&&) = delete;

 public:
  /** destructor */
  virtual ~ControlSimulation() noexcept = default;

  /** Modify the value of theta (for MoreauJeanOSI)
   * \param newTheta the new value of theta */
  void setTheta(double newTheta);

  /** Initialize the ControlSimulation, instantiate all objects
   */
  void initialize();

  /** Add a DynamicalSystem
   * \param ds the DynamicalSystem to integrate
   * \param name of the ds (optional)
   */
  void addDynamicalSystem(std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
                          const std::string& name = "");

  /** Add a Sensor
   * \param sensor the sensor to be added
   * \param h sampling period (or timestep) for the Sensor
   */
  void addSensor(std::shared_ptr<Sensor> sensor, const double h);

  /** Add an Actuator
   * \param actuator the controller to be added
   * \param h sampling period (or timestep) for the Actuator
   */
  void addActuator(std::shared_ptr<Actuator> actuator, const double h);

  /** Add an Observer
   * \param observer the observer to be added
   * \param h sampling period (or timestep) for the Observer
   */
  void addObserver(std::shared_ptr<Observer> observer, const double h);

  /** store the simulation data in a row of the matrix
   * \param indx the current row index
   */
  void storeData(size_t indx);

  /** Return the Simulation
   * \return the simulation for the main simulation
   */
  inline std::shared_ptr<siconos::simulation::Simulation> simulation() const {
    return _processSimulation;
  };

  /** Return the OneStepIntegrator
   * \return the Integrator
   */
  inline std::shared_ptr<siconos::integrators::OneStepIntegrator> integrator() const {
    return _processIntegrator;
  };

  /** \return the siconos::modeling::NonSmoothDynamicalSystem
   */
  inline std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> model() const {
    return _nsds;
  }

  /** Return the data matrix
   * \return the data matrix
   */
  inline std::shared_ptr<siconos::algebra::SiconosMatrix> data() const { return _dataM; }

  /** get the legend for the matrix
   * \return legend as string of space seperated values
   */
  inline std::string dataLegend() const { return _dataLegend; }

  /**
     \return the elapsed time computing in seconds
  */
  inline double elapsedTime() const { return _elapsedTime; }

  /** Return the ControlManager
   * \return the ControlManager
   */
  inline std::shared_ptr<ControlManager> CM() const { return _CM; };

  /** Set the value of _saveOnlyMainSimulation
   * \param v a boolean
   */
  inline void setSaveOnlyMainSimulation(bool v) { _saveOnlyMainSimulation = v; };

  /** Set the simulation to be silent, e.g. do not show any progress bar
   * \param s is true is silent, else display progress bar */
  void silent(bool s = true) { _silent = s; };

  /** Run the simulation */
  virtual void run() = 0;
};
}  // namespace siconos::control

#endif  // CONTROLDYNAMICALSYSTEM_H
