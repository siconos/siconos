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

/*! \file ControlSimulation.hpp
  \brief Abstract class - General interface for all Control Dynamical Systems.
*/

#ifndef CONTROLDYNAMICALSYSTEM_H
#define CONTROLDYNAMICALSYSTEM_H

#include "SiconosPointers.hpp"
#include "SiconosAlgebraTypeDef.hpp"
#include "ControlTypeDef.hpp"
#include "SiconosControlFwd.hpp"
#include "SiconosFwd.hpp"

#include <string>

class ControlSimulation
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(ControlSimulation);

protected:
  /** default constructor */
  ControlSimulation() {};

  /** Constructor with the minimal set of data
   * \param t0 the starting time \f$t_0\f$
   * \param T the end time T
   * \param h the simulation time step
   * */
  ControlSimulation(double t0, double T, double h);

  /** Starting time */
  double _t0;

  /** End time */
  double _T;

  /** Simulation step */
  double _h;

  /** Theta for MoreauJeanOSI */
  double _theta;

  /** Time spent computing */
  double _elapsedTime;

  /** rough estimation of the number of points to save */
  unsigned _N;

  /** Dimension of the state space */
  unsigned _nDim;

  /** Save only the data in the main Simulation*/
  bool _saveOnlyMainSimulation;

  /** If true, do not show progress of the simulation */
  bool _silent;

  /** Matrix for saving result */
  SP::SimpleMatrix _dataM;

  /** Legend for the columns in the matrix _dataM*/
  std::string _dataLegend;

  /** Model */
  SP::Model _model;

  /** TimeDiscretisation for the simulation*/
  SP::TimeDiscretisation _processTD;

  /** The Simulation */
  SP::Simulation _processSimulation;

  /** The integrator */
  SP::OneStepIntegrator _processIntegrator;

  /** the ControlManager */
  SP::ControlManager _CM;

  /** DynamicalSystemsGraph (for convenience)*/
  SP::DynamicalSystemsGraph _DSG0;

  /** InteractionsGraph (for convenience)*/
  SP::InteractionsGraph _IG0;

public:

  /** destructor */
  virtual ~ControlSimulation() {};

  /** Modify the value of theta (for MoreauJeanOSI)
   * \param newTheta the new value of theta */
  void setTheta(unsigned int newTheta);

  /** Initialize the ControlSimulation, instantiate all objects
   */
  void initialize();

  /** Add a DynamicalSystem
   * \param ds the DynamicalSystem to integrate
   * \param name of the ds (optional)
   */
  void addDynamicalSystem(SP::DynamicalSystem ds, const std::string& name = "");

  /** Add a Sensor
   * \param sensor the sensor to be added
   * \param h sampling period (or timestep) for the Sensor
   */
  void addSensor(SP::Sensor sensor, const double h);

  /** Add an Actuator
   * \param actuator the controller to be added
   * \param h sampling period (or timestep) for the Actuator
   */
  void addActuator(SP::Actuator actuator, const double h);

  /** Add an Observer
   * \param observer the observer to be added
   * \param h sampling period (or timestep) for the Observer
   */
  void addObserver(SP::Observer observer, const double h);

  /** store the simulation data in a row of the matrix
   * \param indx the current row index
   */
  void storeData(unsigned indx);

  /** Return the Simulation
   * \return the simulation for the main simulation
  */
  inline SP::Simulation simulation() const
  {
    return _processSimulation;
  };

  /** Return the OneStepIntegrator
   * \return the Integrator
  */
  inline SP::OneStepIntegrator integrator() const
  {
    return _processIntegrator;
  };

  /** Return the Model
   * \return the Model 
  */
  inline SP::Model model() const
  {
    return _model;
  }

  /** Return the data matrix
   * \return the data matrix
   */
  inline SP::SimpleMatrix data() const
  {
    return _dataM;
  }

  /** get the legend for the matrix
   * \return legend as string of space seperated values
   */
  inline std::string dataLegend() const
  {
    return _dataLegend;
  }

  /** Return the elapsed time computing
   * \return the elapsed time computing
   */
  inline double elapsedTime() const
  {
    return _elapsedTime;
  }

  /** Return the ControlManager
   * \return the ControlManager
   */
  inline SP::ControlManager CM() const
  {
    return _CM;
  };

  /** Set the value of _saveOnlyMainSimulation
   * \param v a boolean
   */
  inline void setSaveOnlyMainSimulation(bool v)
  {
    _saveOnlyMainSimulation = v;
  };

  /** Set the simulation to be silent, e.g. do not show any progress bar
   * \param s is true is silent, else display progress bar */
  void silent(bool s = true) { _silent = s; };

  /** Run the simulation */
  virtual void run() = 0;
};

#endif // CONTROLDYNAMICALSYSTEM_H
