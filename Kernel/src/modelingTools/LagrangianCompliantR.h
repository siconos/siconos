/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
/*! \file LagrangianCompliantR.h

*/
#ifndef LagrangianCompliantR_H
#define LagrangianCompliantR_H

#include "LagrangianR.h"

class DynamicalSystem;

/** Lagrangian (Non Linear) Compliant Relation: Scleronomous, Non-Holonomic (function of lambda).
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date Apr 27, 2004
 *
 * \f[
 * Y[0] = y = h(q,\lambda(t),z)
 * \f]
 *
 * \f[
 * Y[1] = \dot y = G0(q,\lambda(t),z)\dot q + G1((q,\lambda(t),z)\dot\lambda(t)
 * \f]
 *
 * \f[
 * p = G0^t(q,\lambda(t),z)\lambda(t)
 * \f]
 *
 * with
 * \f[
 * G0(q,\lambda(t),z) = \nabla_q h(q,\lambda(t),z)
 * \f]
 * \f[
 * G1(q,\lambda(t),z) = \nabla_{\lambda}h(q,\lambda(t),z)
 * \f]
 *
 * h, G0 and G1 are connected to user-defined functions.
 *
 */
class LagrangianCompliantR : public LagrangianR
{

protected:
  /** LagrangianR plug-in to compute h(q,lambda,z)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param sizeY : size of vector y (ie of lambda and of the interaction)
   * @param lambda : pointer to lambda of the interaction
   * @param[in,out] y : pointer to the first element of y
   * @param sizeZ : size of vector z.
   * @param[in,out] z : a vector of user-defined parameters
   */
  FPtr2 hPtr;

  /** LagrangianR plug-in to compute G0(q,lambda,z)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param sizeY : size of vector y (ie of lambda and of the interaction)
   * @param lambda : pointer to lambda of the interaction
   * @param[in,out] G0 : pointer to the first element of G0
   * @param sizeZ : size of vector z.
   * @param[in,out] z : a vector of user-defined parameters
   */
  FPtr2 G0Ptr;

  /** LagrangianR plug-in to compute G1(q,lambda,z)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param sizeY : size of vector y (ie of lambda and of the interaction)
   * @param lambda : pointer to lambda of the interaction
   * @param[in,out] G1 : pointer to the first element of G1
   * @param sizeZ : size of vector z.
   * @param[in,out] z : a vector of user-defined parameters
   */
  FPtr2 G1Ptr;

public:

  /** default constructor
   */
  LagrangianCompliantR();

  /** constructor from xml file
   *  \param relationXML
   */
  LagrangianCompliantR(RelationXML*);

  /** constructor from a set of data
   *  \param string : the name of the plugin to computeH
   *  \param vector<string> : a list of names for the plugin to computeG (depends on relation type)
   */
  LagrangianCompliantR(const std::string&, const std::vector<std::string>&);

  /** destructor
   */
  ~LagrangianCompliantR();

  /** initialize G matrices or components specific to derived classes.
   */
  void initComponents();

  /** to set a specified function to compute function h(q,...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeHFunction(const std::string& , const std::string&);

  /** to set a specified function to compute G(q, ...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \param unsigned int: the index of G that must be computed (see introduction of this class for details on indexes)
   *  \exception SiconosSharedLibraryException
   */
  void setComputeGFunction(const std::string& , const std::string&  , unsigned int  = 0);

  /** to compute y = h(q,v,t) using plug-in mechanism
   * \param: double, current time
   */
  void computeH(double);

  /** to compute G using plug-in mechanism. Index shows which G is to be computed
   * \param: double, current time
   * \param: unsigned int
   */
  void computeG(double, unsigned int = 0);

  /** to compute output
   *  \param Interaction : the interaction that owns y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(double, unsigned int = 0);

  /** to compute p
   *  \param Interaction : the interaction that owns lambda
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  void computeInput(double, unsigned int);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation* : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianCompliantR* convert(Relation *r);
};

#endif
