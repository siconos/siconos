/* Siconos-Kernel, Copyright INRIA 2005-2013.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */


#ifndef MatrixIntegrator_H
#define MatrixIntegrator_H

/* ! \file MatrixIntegrator.hpp
 * \brief Class to integrate a Matrix ODE \f$\dot{X} = AX + E\f$
 */

#include "SiconosSerialization.hpp"
#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosPointers.hpp"
#include "OneStepIntegratorTypes.hpp"


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

  /** Model for integration */
  SP::Model _model;

  /** TimeDiscretisation for integration */
  SP::TimeDiscretisation _TD;

  /** Simulation (of EventDriven type) */
  SP::EventDriven _sim;

  /** OneStepIntegrator of type Lsodar */
  SP::Lsodar _OSI;

  /** */
  void commonInit(const DynamicalSystem& ds, const Model& m);

  /** Default constructor */
  MatrixIntegrator() {};

public:

  /** Constructor to compute \f$\int exp(A\tau)E\amthrm{d}\tau\f$
   * \param ds the DynamicalSystem
   * \param m the original Model
   * \param E a matrix
   */
  MatrixIntegrator(const DynamicalSystem& ds, const Model& m, SP::SiconosMatrix E);

  /** Constructor to compute \f$\int exp(A\tau)E(\tau)\mathrm{d}\tau\f$
   * \param ds the DynamicalSystem
   * \param m the original Model
   * \param plugin the plugin to compute \f$E(t)\f$
   * \param p the number of column in E
   */
  MatrixIntegrator(const DynamicalSystem& ds, const Model& m, SP::PluggedObject plugin, const unsigned int p);

  /** Constructor to compute \f$\int exp(A\tau)\mathrm{d}\tau\f$
   * \param ds the DynamicalSystem
   * \param m the original Model
   */
  MatrixIntegrator(const DynamicalSystem& ds, const Model& m);

  /** Computes the next value of _mat */
  void integrate();

  /** Get the value of _mat, solution of the ODE
   * \return a reference to _mat */
  inline const SiconosMatrix& mat() const { return *_mat; }

  /** Check whether the solution of the ODE is time-invariant*/
  inline bool isConst() { return _isConst; }

};

#endif
