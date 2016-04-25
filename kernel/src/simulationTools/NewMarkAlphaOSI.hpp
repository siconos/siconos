/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
/*! \file
  NewMark Alpha Scheme Time-Integrator for Dynamical Systems
*/
#ifndef NEWMARKALPHAOSI_H
#define NEWMARKALPHAOSI_H

#include "OneStepIntegrator.hpp"

/**  NewMarkAlpha Scheme Time-Integrator for Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * NewMarkAlphaOSI is used to solve constrained dynamical systems represented by index-3 DAE
 *
 * NewMarkAlphaOSI is instantiated with values of beta, gamma, alpha_m, alpha_f and the list of concerned
 * dynamical systems. Each DynamicalSystem is associated to a SiconosMatrix named "W"
 *
 * W matrices are initialized and computed in initW and computeW.
 */
class NewMarkAlphaOSI : public OneStepIntegrator
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewMarkAlphaOSI);
  /** Parameters of the numerical scheme:  beta, gamma, alpha_m, alpha_f */
  double _beta, _gamma, _alpha_m, _alpha_f;

  /** Stl map that associates a W matrix to each DynamicalSystem of the OSI */
  MapOfDSMatrices WMap;

  /** Order of the polynomial for dense output*/
  unsigned int _orderDenseOutput;

  /** Indicator whether or not constraints at the velocity level are handled
   * _IsVelocityLevel = true: constraints at the velocity level are handled
   * _IsVelocityLevel = false: constraints at the position are handled
   */
  bool _IsVelocityLevel;
  /**
   * Default constructor
  */
  NewMarkAlphaOSI() {};

public:
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
  virtual ~NewMarkAlphaOSI() {};


  // --- GETTERS/SETTERS ---

  /** set value to the parameter beta
   * \param beta value of beta
   */
  inline void setBeta(double beta)
  {
    _beta = beta;
  };

  /** set value to the parameter gamma
   * \param value_gamma double : value of gamma
   */
  inline void setGamma(double value_gamma)
  {
    _gamma = value_gamma;
  };

  /** set value to the parameter alpha_m
   * \param value_alpha_m double : value of alpha_m
   */

  inline void setAlpha_m(double value_alpha_m)
  {
    _alpha_m = value_alpha_m;
  };

  /** set value to the parameter alpha_f
   * \param value_alpha_f double : value of alpha_f
   */

  inline void setAlpha_f(double value_alpha_f)
  {
    _alpha_f = value_alpha_f;
  };

  /** set values to the parameters beta, gamma, alpha_f, alpha_m from the value of rho_infty
   * \param rho_infty double : value of rho_infty
   */

  inline void setParametersFromRho_infty(double rho_infty)
  {
    _alpha_m = (2 * rho_infty - 1) / (rho_infty + 1);
    _alpha_f = rho_infty / (rho_infty + 1);
    _gamma = 0.5 + _alpha_f - _alpha_m;
    _beta = 0.25 * std::pow((_gamma + 0.5), 2);
  };

  /** get value of beta
   * \return double
   */
  inline double getBeta()
  {
    return _beta;
  };

  /** get value of gamma
   * \return double
   */
  inline double getGamma()
  {
    return _gamma;
  };

  /** get value of alpha_m
   * \return double
   */
  inline double getAlpha_m()
  {
    return _alpha_m;
  };

  /** get value of alpha_f
   * \return double
   */

  inline double getAlpha_f()
  {
    return _alpha_f;
  };
  /** get the order of the polynomial for dense output
   * \return unsigned int
   */
  inline unsigned int getOrderDenseOutput()
  {
    return _orderDenseOutput;
  }

  /** set the flag _IsVelocityLevel
   * \param flag bool
   */
  inline void setFlagVelocityLevel(bool flag)
  {
    _IsVelocityLevel = flag;
  }

  /** get the flag _IsVelocityLevel
   * \return bool
   */
  inline bool getFlagVelocityLevel()
  {
    return _IsVelocityLevel;
  }

  /** get matrix W
   * \param ds SP::DynamicalSystem DynamicalSystem concerned
   * \return  SimpleMatrix
   */
  const SimpleMatrix getW(SP::DynamicalSystem ds);

  /** get pointer to the maxtrix W
   * \param ds SP::DynamicalSystem DynamicalSystem concerned
   * \return  SP::SimpleMatrix
   */
  SP::SimpleMatrix W(SP::DynamicalSystem ds);

  /** initialize WMap[ds] matrix
    *  \param ds a pointer to DynamicalSystem
    */
  void initW(SP::DynamicalSystem ds );

  /** compute WMap[ds] matrix
   *  \param ds a pointer to DynamicalSystem
   */
  void computeW(SP::DynamicalSystem ds);

  /** compute the residual of dynamical equation
   *\return double: maximum residu over all DSs
   */
  double computeResidu();

  /** compute the free state of the discretized dynamical system */
  void computeFreeState();

  /** integrates the Interaction linked to this integrator, without taking non-smooth effects into account
   * \param vertex_inter of the interaction graph
   * \param osnsp pointer to OneStepNSProblem
   */
  virtual void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter,
                                 OneStepNSProblem* osnsp);
  
  /** initialize */
  void initialize(Model& m);

  /** prepare for Newton Iteration 
   * \param time
   */
  void prepareNewtonIteration(double time);

  /** predict first values for the Newton iteration */
  void prediction();

  /** correct state of all levels of Dynamical Systems after each Newton iteration
   */
  void correction();

  /** integrate the system, between tinit and tend with possible stop at tout
   *  \param tinit double: tinit, initial time
   *  \param tend double: tend, end time
   *  \param tout double: tout, real end time
   *  \param flag useless for NewMarkAlphaOSI
   */
  void integrate(double& tinit, double& tend, double& tout, int& flag);

  /** updates the state of the Dynamical Systems
   *  \param level the level of interest for the dynamics: not used at the time
   */
  void updateState(const unsigned int level);

  /** Compute coefficients of the polynomial of the dense output for a given DS
   *  \param ds SP::DynamicalSystem, ds concerned
   */
  void computeCoefsDenseOutput(SP::DynamicalSystem ds);

  /** prepare for Event localization*/
  void prepareEventLocalization();

  /** Generate dense output for all Dynamical Systems belonging to OSI
   * \param time at which we want to generate the dense output
   */
  void DenseOutputallDSs(double time);

  /** Displays the data of the NewMarkAlpha's integrator
   */
  void display();






  ACCEPT_STD_VISITORS();

};

#endif // NEWMARKALPHAOSI_H
