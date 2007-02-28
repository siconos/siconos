/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
/*! \file LagrangianRheonomousR.h

*/
#ifndef LagrangianRheonomousR_H
#define LagrangianRheonomousR_H

#include "LagrangianR.h"

class DynamicalSystem;

/** Lagrangian (Non Linear) Relation, Rheonomous and Holonomic.
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date February 28, 2007
 *
 *  This class provides tools to describe non linear relation of the type:
 *
 * \f[
 * Y[0] = y = h(q,t,z)
 * \f]
 *
 * \f[
 * Y[1] = \dot y = G0(q,t,z)\dot q + \frac{\partial h}{\partial t}(q,t,z)
 * \f]
 *
 * \f[
 * p = G0^t(q,t,z)\lambda
 * \f]
 *
 * with
 * \f[
 * G0(q,\lambda,z) = \nabla_q h(q,\lambda,z)
 * \f]
 *
 *
 * h, G0 and hdot=\f$ \frac{\partial h}{\partial t}(q,t,z) \f$ are connected to user-defined functions.
 *
 */
class LagrangianRheonomousR : public LagrangianR
{

protected:

  /** derivative of h according to time */
  SiconosVector * hDot;

  /** LagrangianRheonomousR plug-in to compute h(q,t,z)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param time : current time
   * @param sizeY : size of vector y (ie of the intercation)
   * @param[in,out] y : pointer to the first element of y
   * @param sizeZ : size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  FPtr4 hPtr;

  /** LagrangianRheonomousR plug-in to compute hDot(q,t,z)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param time : current time
   * @param sizeY : size of vector hDot (ie of the intercation)
   * @param[in,out] pointer to the first element of hDot
   * @param sizeZ : size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  FPtr4 hDotPtr;

  /** LagrangianRheonomousR plug-in to compute G0(q,t,z), gradient of h accoring to q
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param time : current time
   * @param sizeY : size of vector y (ie of the intercation)
   * @param[in,out] G0 : pointer to the first element of G0
   * @param sizeZ : size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  FPtr4 G0Ptr;

public:

  /** default constructor
   */
  LagrangianRheonomousR();

  /** constructor from xml file
   *  \param relationXML
   *  \exception RuntimeException
   */
  LagrangianRheonomousR(RelationXML*);

  /** constructor from a set of data
   *  \param string : the name of the plugin to compute h
   *  \param string : the name of the plugin to compute hDot
   *  \param string : the name of the plugin to computeG0
   */
  LagrangianRheonomousR(const std::string&, const std::string&, const std::string&);

  /** destructor
   */
  ~LagrangianRheonomousR();

  /** initialize G matrices or components specific to derived classes.
   */
  void initComponents();

  /** get vector hDot
   *  \return a SimpleVector.
   */
  inline const SimpleVector getHDot() const
  {
    return *hDot;
  }

  /** get a pointer on vector hDot
   *  \return a pointer on a SiconosVector.
   */
  inline SiconosVector* getHDotPtr() const
  {
    return hDot;
  }

  /** set the value of hDot to newValue (copy)
   *  \param a SiconosVector.
   */
  void setHDot(const SiconosVector&);

  /** set hDot to pointer newPtr (pointer link)
   *  \param  SiconosVector* newPtr
   */
  void setHDotPtr(SiconosVector *newPtr);

  /** to set a specified function to compute function h(q,...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeHFunction(const std::string& , const std::string&);

  /** to set a specified function to compute function hDot
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeHDotFunction(const std::string& , const std::string&);

  /** to set a specified function to compute G(q, ...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \param unsigned int: the index of G that must be computed (see introduction of this class for details on indexes)
   */
  void setComputeGFunction(const std::string& , const std::string&  , unsigned int  = 0);

  /** to compute y = h(q,v,t) using plug-in mechanism
   * \param: double, current time
   */
  void computeH(double);

  /** to compute hDot using plug-in mechanism
   * \param: double, current time
   */
  void computeHDot(double);

  /** to compute G using plug-in mechanism. Index shows which G is to be computed
   * \param: double, current time
   * \param: unsigned int
   */
  void computeG(double, unsigned int = 0);

  /** to compute y = h(q,v,t) using plug-in mechanism, with free DS states.
   * \param: double, current time
   */
  void computeHFree(double);

  /** to compute hDot using plug-in mechanism, with free DS states.
   * \param: double, current time
   */
  void computeHDotFree(double);

  /** to compute G using plug-in mechanism, with free DS states.. Index shows which G is to be computed
   * \param: double, current time
   * \param: unsigned int
   */
  void computeGFree(double, unsigned int = 0);

  /** to compute output
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double, unsigned int = 0);

  /** to compute y for the free state
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeFreeOutput(double, unsigned int = 0);

  /** to compute p
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  void computeInput(double, unsigned int);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation* : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianRheonomousR* convert(Relation *r);
};

#endif
