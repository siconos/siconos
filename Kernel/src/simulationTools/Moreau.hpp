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
  Moreau Time-Integrator for Dynamical Systems
*/

#ifndef MOREAU_H
#define MOREAU_H

#include "OneStepIntegrator.hpp"
#include "SimpleMatrix.hpp"

class Simulation;
class SiconosMatrix;

const unsigned int MOREAUSTEPSINMEMORY = 1;

/**  Moreau Time-Integrator for Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * See User's guide, \ref docSimuMoreauTS for details.
 *
 * Moreau class is used to define some time-integrators methods for a
 * list of dynamical systems.

 * A Moreau instance is defined by the value of theta and the list of
 * concerned dynamical systems.  Each DynamicalSystem is associated to
 * a SiconosMatrix, named "W"
 *
 * W matrices are initialized and computed in initW and
 * computeW. Depending on the DS type, they may depend on time and DS
 * state (x).
 *
 * For first order systems, the implementation uses _r for storing the
 * the input due to the nonsmooth law. This Moreau scheme assumes that the
 * relative degree is zero or one and one level for _r is sufficient
 *
 * For Lagrangian systems, the implementation uses _p[1] for storing the
 * discrete impulse.
 *
 * Main functions:
 *
 * - computeFreeState(): computes xfree (or vfree), dynamical systems
 *   state without taking non-smooth part into account \n
 *
 * - updateState(): computes x (q,v), the complete dynamical systems
 *    states.
 *
 */
class Moreau : public OneStepIntegrator
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Moreau);


  /** Stl map that associates a W Moreau matrix to each DynamicalSystem of the OSI */
  MapOfDSMatrices WMap;

  /** Stl map that associates the columns of  W Moreau matrix to each DynamicalSystem of the OSI if it has some boundary conditions */
  MapOfDSMatrices _WBoundaryConditionsMap;

  /** Stl map that associates a theta parameter for the integration
  *  scheme to each DynamicalSystem of the OSI */
  double _theta;

  /** A gamma parameter for the integration scheme to each DynamicalSystem of the OSI
   * This parameter is used to apply a theta-method to the input $r$
   */
  double _gamma;

  /** a boolean to known is the parameter must be used or not
   */
  bool _useGamma;

  /** a boolean to known is the parameter must be used or not
   */
  bool _useGammaForRelation;

  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;



  /** Default constructor
   */
  Moreau() {};

public:

  /** constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param DynamicalSystemsSet: set of all DS in the NSDS
   */
  Moreau(SP::OneStepIntegratorXML, SP::DynamicalSystemsSet);

  /** constructor from a minimum set of data: one DS and its theta
   *  \param ds SP::DynamicalSystem : the DynamicalSystem linked to the OneStepIntegrator
   *  \param theta value of the parameter
   */
  Moreau(SP::DynamicalSystem, double);

  /** constructor from a minimum set of data
   *  \param DynamicalSystemsSet : the list of DynamicalSystems to be integrated
   *  \param theta value for all these DS.
   */
  Moreau(DynamicalSystemsSet&, double);

  /** constructor from theta value only
   *  \param theta value for all these DS.
   */
  Moreau(double);

  /** constructor from a minimum set of data: one DS and its theta
   *  \param SP::DynamicalSystem : the DynamicalSystem linked to the OneStepIntegrator
   *  \param Theta value
   *  \param gamma value
   */
  Moreau(SP::DynamicalSystem, double, double);

  /** constructor from a minimum set of data
   *  \param DynamicalSystemsSet : the list of DynamicalSystems to be integrated
   *  \param gamma value for all these DS.
   *  \param theta value for all these DS.
   */
  Moreau(DynamicalSystemsSet&, double, double);

  /** constructor from theta value only
   *  \param theta value for all these DS.
   *  \param gamma value for all these DS.
   */
  Moreau(double, double);

  /** destructor
   */
  virtual ~Moreau() {};

  // --- GETTERS/SETTERS ---

  /** get the value of W corresponding to DynamicalSystem ds
   * \param a pointer to DynamicalSystem, optional, default =
   * NULL. get W[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getW(SP::DynamicalSystem = SP::DynamicalSystem());

  /** get W corresponding to DynamicalSystem ds
   * \param a pointer to DynamicalSystem, optional, default =
   * NULL. get W[0] in that case
   * \return pointer to a SiconosMatrix
   */
  SP::SimpleMatrix W(SP::DynamicalSystem ds);

  /** set the value of W[ds] to newValue
   * \param SiconosMatrix newValue
   * \param a pointer to DynamicalSystem,
   */
  void setW(const SiconosMatrix&, SP::DynamicalSystem);

  /** set W[ds] to pointer newPtr
   * \param SP::SiconosMatrix  newPtr
   * \param a pointer to DynamicalSystem
   */
  void setWPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem);

  // -- WBoundaryConditions --

  /** get WBoundaryConditions map
   *  \return a MapOfDSMatrices
   */
  inline MapOfDSMatrices getWBoundaryConditionsMap() const
  {
    return _WBoundaryConditionsMap;
  };

  /** get the value of WBoundaryConditions corresponding to DynamicalSystem ds
   * \param a pointer to DynamicalSystem, optional, default =
   * NULL. get WBoundaryConditions[0] in that case
   *  \return SimpleMatrix
   */
  const SimpleMatrix getWBoundaryConditions(SP::DynamicalSystem = SP::DynamicalSystem());

  /** get WBoundaryConditions corresponding to DynamicalSystem ds
   * \param a pointer to DynamicalSystem, optional, default =
   * NULL. get WBoundaryConditions[0] in that case
   * \return pointer to a SiconosMatrix
   */
  SP::SiconosMatrix WBoundaryConditions(SP::DynamicalSystem ds);

  // -- theta --

  /** get theta
   *  \return a double
   */
  inline double theta()
  {
    return _theta;
  };

  /** set the value of theta
   *  \param a double
   */
  inline void setTheta(double newTheta)
  {
    _theta = newTheta;
  };

  // -- gamma --

  /** get gamma
   *  \return a double
   */
  inline double gamma()
  {
    return _gamma;
  };

  /** set the value of gamma
   *  \param a double
   */
  inline void setGamma(double newGamma)
  {
    _gamma = newGamma;
    _useGamma = true;
  };

  // -- useGamma --

  /** get bool useGamma
   *  \return a bool
   */
  inline bool useGamma()
  {
    return _useGamma;
  };

  /** set the boolean to indicate that we use gamma
   *  \param a bool
   */
  inline void setUseGamma(bool newUseGamma)
  {
    _useGamma = newUseGamma;
  };

  /** get bool gammaForRelation for the relation
   *  \return a
   */
  inline bool useGammaForRelation()
  {
    return _useGammaForRelation;
  };


  /** set the boolean to indicate that we use gamma for the relation
   *  \param a bool
   */
  inline void setUseGammaForRelation(bool newUseGammaForRelation)
  {
    _useGammaForRelation = newUseGammaForRelation;
    if (_useGammaForRelation) _useGamma = false;
  };


  // --- OTHER FUNCTIONS ---

  /** initialization of the Moreau integrator; for linear time
      invariant systems, we compute time invariant operator (example :
      W)
   */
  virtual void initialize();

  /** init WMap[ds] Moreau matrix at time t
   *  \param the time (double)
   *  \param a pointer to DynamicalSystem
   */
  void initW(double, SP::DynamicalSystem);

  /** compute WMap[ds] Moreau matrix at time t
   *  \param the time (double)
   *  \param a pointer to DynamicalSystem
   */
  void computeW(double, SP::DynamicalSystem);

  /** compute WBoundaryConditionsMap[ds] Moreau matrix at time t
   *  \param the time (double)
   *  \param a pointer to DynamicalSystem
   */
  void computeWBoundaryConditions(SP::DynamicalSystem);

  /** init WBoundaryConditionsMap[ds] Moreau
   *  \param a pointer to DynamicalSystem
   */
  void initWBoundaryConditions(SP::DynamicalSystem);

  /** return the maximum of all norms for the "Moreau-discretized" residus of DS
      \return a double
   */
  double computeResidu();

  /** Perform the integration of the dynamical systems linked to this integrator
   *  without taking into account the nonsmooth input (_r or _p)
   */
  virtual void computeFreeState();

  /** Compute the Output (y) which corresponds to the free state (state without
   * taking into account the nonsmooth input) plus the possible contribution of
   * the nslaw
   */
  virtual void computeFreeOutput(SP::Interaction inter, OneStepNSProblem * osnsp);

  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   */
  virtual bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   */
  virtual bool removeInteractionInIndexSet(SP::Interaction inter, unsigned int i);


  /**
   */
  void prepareNewtonIteration(double time);


  /** integrate the system, between tinit and tend (->iout=true), with possible stop at tout (->iout=false)
   *  \param double: tinit, initial time
   *  \param double: tend, end time
   *  \param double: tout, real end time
   *  \param int: useless flag (for Moreau, used in Lsodar)
   */
  void integrate(double&, double&, double&, int&);

  /** updates the state of the Dynamical Systems
   *  \param level the level of interest for the dynamics: not used at the time
   */
  virtual void updateState(const unsigned int level);

  /** copy the matrix W of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveWToXML();

  /** Displays the data of the Moreau's integrator
   */
  void display();

  /** insert a dynamical system in this Integrator
   *  \param a SP::DynamicalSystem
   */
  void insertDynamicalSystem(SP::DynamicalSystem ds);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, 0 otherwise
   */
  static Moreau* convert(OneStepIntegrator* osi);

  /** visitors hook
  */
  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(Moreau)

#endif // MOREAU_H
