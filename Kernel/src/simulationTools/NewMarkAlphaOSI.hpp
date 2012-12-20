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
using namespace std;
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
  double beta, gamma, alpha_m, alpha_f;

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

  /** constructor with one DS and parameters beta, gamma, alpha_m, alpha_f
   *  \param ds the DynamicalSystem linked to the OneStepIntegrator
   *  \param beta value of beta
   *  \param gamma double : value of gamma
   *  \param alpha_m double : value of alpha_m
   *  \param alpha_f double : value of alpha_f
   */
  NewMarkAlphaOSI(SP::DynamicalSystem, double, double, double, double, bool);

  /** constructor with one DS and the parameter rho_infty from which values of beta, gamma, alpha_m, alpha_f can be deduced
   * \param ds SP::DynamicalSystem
   * \param rho_infty double : value of rho_infty
   */

  NewMarkAlphaOSI(SP::DynamicalSystem, double, bool);

  /** constructor with a list of concerned Dynamical Systems and parameters beta, gamma, alpha_m, alpha_f
  * \param listOfDS DynamicalSystemsSet : list of Dynamical Systems to be integrated
  * \param beta double
  * \param gamma double
  * \param alpha_m double
  * \param alpha_f double
  */

  NewMarkAlphaOSI(DynamicalSystemsSet&, double, double, double, double, bool);

  /** constructor with a list of concerned Dynamical Systems and the parameter rho_infty
   * \param listOfDS DynamicalSystemsSet
   * \param rho_infty double
   */

  NewMarkAlphaOSI(DynamicalSystemsSet&, double, bool);

  /** constructor with only parameters beta, gamma, alpha_m, alpha_f
  * \param beta double
  * \param gamma double
  * \param alpha_m double
  * \param alpha_f double
  */

  NewMarkAlphaOSI(double, double, double, double, bool);

  /** constructor with only the parameter rho_infty
  * \param rho_infty double
  */

  NewMarkAlphaOSI(double, bool);

  /** destructor
   */
  virtual ~NewMarkAlphaOSI() {};


  // --- GETTERS/SETTERS ---

  /** set value to the parameter beta
   * \param double : value of beta
   */

  inline void setBeta(double value_beta)
  {
    beta = value_beta;
  };

  /** set value to the parameter gamma
   * \param double : value of gamma
   */

  inline void setGamma(double value_gamma)
  {
    gamma = value_gamma;
  };

  /** set value to the parameter alpha_m
   * \param double : value of alpha_m
   */

  inline void setAlpha_m(double value_alpha_m)
  {
    alpha_m = value_alpha_m;
  };

  /** set value to the parameter alpha_f
   * \param double : value of alpha_f
   */

  inline void setAlpha_f(double value_alpha_f)
  {
    alpha_f = value_alpha_f;
  };

  /** set values to the parameters beta, gamma, alpha_f, alpha_m from the value of rho_infty
   * \param double : value of rho_infty
   */

  inline void setParametersFromRho_infty(double _rho_infty)
  {
    alpha_m = (2 * _rho_infty - 1) / (_rho_infty + 1);
    alpha_f = _rho_infty / (_rho_infty + 1);
    gamma = 0.5 + alpha_f - alpha_m;
    beta = 0.25 * std::pow((gamma + 0.5), 2);
  };

  /** get value of beta */

  inline double getBeta()
  {
    return beta;
  };

  /** get value of gamma*/

  inline double getGamma()
  {
    return gamma;
  };

  /** get value of alpha_m*/

  inline double getAlpha_m()
  {
    return alpha_m;
  };

  /** get value of alpha_f */

  inline double getAlpha_f()
  {
    return alpha_f;
  };
  /** get the order of the polynomial for dense output */
  inline unsigned int getOrderDenseOutput()
  {
    return _orderDenseOutput;
  }

  /** set the flag _IsVelocityLevel
   * \param bool
   */
  inline void setFlagVelocityLevel(bool _flag)
  {
    _IsVelocityLevel = _flag;
  }

  /** get the flag _IsVelocityLevel
   * \return bool
   */
  inline bool getFlagVelocityLevel()
  {
    return _IsVelocityLevel;
  }

  /** get matrix W
   *\param SP::DynamicalSystem DynamicalSystem concerned
   */
  const SimpleMatrix getW(SP::DynamicalSystem ds);

  /** get pointer to the maxtrix W
   *\param SP::DynamicalSystem DynamicalSystem concerned
   */

  SP::SimpleMatrix W(SP::DynamicalSystem ds);

  /** insert a dynamical system in this Integrator
   *  \param a SP::DynamicalSystem
   */
  void insertDynamicalSystem(SP::DynamicalSystem ds);


  /** initialize WMap[ds] matrix
    *  \param a pointer to DynamicalSystem
    */
  void initW(SP::DynamicalSystem);

  /** compute WMap[ds] matrix
   *  \param a pointer to DynamicalSystem
   */
  void computeW(SP::DynamicalSystem);

  /** compute the residual of dynamical equation
   *\return double: maximum residu over all DSs
   */
  double computeResidu();

  /** compute the free state of the discretized dynamical system */
  void computeFreeState();

  /** compute free output (without interaction force) of the concerned interaction
  /param SP::Interaction: pointer to the concerned interaction
  /param OneStepNSProblem: method to solve non-smooth problem
  */
  void computeFreeOutput(SP::Interaction, OneStepNSProblem *);

  /** initialize */
  void initialize();

  /** prepare for Newton Iteration */
  void prepareNewtonIteration(double);

  /** predict first values for the Newton iteration */
  void prediction();

  /** correct state of all levels of Dynamical Systems after each Newton iteration
   */
  void correction();

  /** integrate the system, between tinit and tend with possible stop at tout
   *  \param double: tinit, initial time
   *  \param double: tend, end time
   *  \param double: tout, real end time
   *  \param int: a flag, useless for NewMarkAlphaOSI
   */
  void integrate(double&, double&, double&, int&);

  /** updates the state of the Dynamical Systems
   *  \param level the level of interest for the dynamics: not used at the time
   */
  void updateState(const unsigned int level);

  /** Compute coefficients of the polynomial of the dense output for a given DS
   *  \param SP::DynamicalSystem, ds concerned
   */
  void computeCoefsDenseOutput(SP::DynamicalSystem);

  /** prepare for Event localization*/
  void prepareEventLocalization();

  /** Generate dense output for all Dynamical Systems belonging to OSI
   * \param double time at which we want to generate the dense output
   */
  void DenseOutputallDSs(double);

  /** Displays the data of the NewMarkAlpha's integrator
   */
  void display();






  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(NewMarkAlphaOSI);
#endif // NEWMARKALPHAOSI_H
