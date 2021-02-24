/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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


#ifndef MatrixIntegrator_H
#define MatrixIntegrator_H

/* ! \file MatrixIntegrator.hpp
 * \brief Class to integrate a Matrix ODE \f$\dot{X} = AX + E\f$
 */

#include "SiconosSerialization.hpp"
#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosPointers.hpp"

#include "SiconosFwd.hpp"

class MatrixIntegrator
{
private:
  /** serialization hooks */
  ACCEPT_SERIALIZATION(MatrixIntegrator);

protected:

  /** The matrix solution to the ODE */
  SP::SiconosMatrix _mat;

  /**The entry Matrix E */
  SP::SiconosMatrix _E;

  /**The entry Matrix E, in the plugin form */
  SP::PluggedObject _plugin;

  /**Plugin to compute the value of a column of E */
  SP::SubPluggedObject _spo;

  /** flag to indicate where the Matrix _mat is constant */
  bool _isConst;

  /** DynamicalSystem to integrate */
  SP::DynamicalSystem _DS;

  /** NonSmoothDynamicalSystem for integration */
  SP::NonSmoothDynamicalSystem _nsds;

  /** TimeDiscretisation for integration */
  SP::TimeDiscretisation _TD;

  /** Simulation (of EventDriven type) */
  SP::EventDriven _sim;

  /** OneStepIntegrator of type LsodarOSI */
  SP::LsodarOSI _OSI;

  /** */
  void commonInit(const DynamicalSystem& ds, const NonSmoothDynamicalSystem& nsds, const TimeDiscretisation & td);

  /** Default constructor */
  MatrixIntegrator() {};

public:

  /** Constructor to compute \f$\int exp(A\tau)E\amthrm{d}\tau\f$
   * \param ds the DynamicalSystem
   * \param nsds current nonsmooth dynamical system
   * \param td current time discretisation
   * \param E a matrix
   */
  MatrixIntegrator(const DynamicalSystem& ds, const NonSmoothDynamicalSystem& nsds, const TimeDiscretisation &td, const  SP::SiconosMatrix E);

  /** Constructor to compute \f$\int exp(A\tau)E(\tau)\mathrm{d}\tau\f$
   * \param ds the DynamicalSystem
   * \param nsds current nonsmooth dynamical system
   * \param td current time discretisation
   * \param plugin the plugin to compute \f$E(t)\f$
   * \param p the number of column in E
   */
  MatrixIntegrator(const DynamicalSystem& ds, const NonSmoothDynamicalSystem& nsds, const TimeDiscretisation &td, SP::PluggedObject plugin, const unsigned int p);

  /** Constructor to compute \f$\int exp(A\tau)\mathrm{d}\tau\f$
   * \param ds the DynamicalSystem
   * \param nsds current nonsmooth dynamical system
   * \param td current time discretisation
   */
  MatrixIntegrator(const DynamicalSystem& ds,const NonSmoothDynamicalSystem& nsds, const TimeDiscretisation &td);

  /** Computes the next value of _mat */
  void integrate();

  /** Get the value of _mat, solution of the ODE
   * \return a reference to _mat */
  inline const SiconosMatrix& mat() const { return *_mat; }

  /** Check whether the solution of the ODE is time-invariant*/
  inline bool isConst() { return _isConst; }

};

#endif
