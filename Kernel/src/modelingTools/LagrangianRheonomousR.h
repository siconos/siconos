/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

/** Lagrangian (Non Linear) Relation, Rheonomous and Holonomic.
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
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
 * G0(q,t,z) = \nabla_q h(q,t,z)
 * \f]
 *
 *
 * h, G0 and hdot=\f$ \frac{\partial h}{\partial t}(q,t,z) \f$ are connected to user-defined functions.
 *
 */
class LagrangianRheonomousR : public LagrangianR<FPtr4>
{
  typedef PluggedObject<FPtr4, SimpleVector> PVT2;
  typedef boost::shared_ptr<PVT2> SPPVT2;

protected:

  typedef LagrangianR<FPtr4> BaseClass ;

  /** plugged vector used to compute hDot */
  SPPVT2 hDot;

  /** initialize G matrices or components specific to derived classes.
   */
  void initComponents();

  /** default constructor
   */
  LagrangianRheonomousR(): BaseClass(RELATION::RheonomousR)
  {
    JacH.resize(1);
  };

public:

  /** constructor from xml file
   *  \param relationXML
   */
  LagrangianRheonomousR(SP::RelationXML);

  /** constructor from a set of data
   *  \param string : the name of the plugin to compute h
   *  \param string : the name of the plugin to compute hDot
   *  \param string : the name of the plugin  to compute jacobian h according to q
   */
  LagrangianRheonomousR(const std::string&, const std::string&, const std::string&);

  /** destructor
   */
  virtual ~LagrangianRheonomousR() {};

  // -- hDot --

  /** get vector hDot
   *  \return a SimpleVector
   */
  inline const SimpleVector getHDot() const
  {
    return *hDot;
  }

  /** get a pointer on vector hDot
   *  \return a smart pointer on a SiconosVector
   */
  inline SP::SiconosVector getHDotPtr() const
  {
    return hDot;
  }

  /** set the value of hDot to newValue (copy)
   *  \param SiconosVector newValue
   */
  template <class U> void setHDot(const U& newValue)
  {
    setObject<PVT2, SPPVT2, U>(hDot, newValue);
  }

  /** set hDot to pointer newPtr (pointer link)
   *  \param SP::SiconosVector  newPtr
   */
  inline void setHDot(SPPVT2 newPtr)
  {
    hDot = newPtr ;
  }

  /** To get the name of hDot plugin
   *  \return a string
   */
  inline const std::string getHDotName() const
  {
    return hDot->getPluginName();
  }

  /** true if hDot is plugged
   *  \return a bool
   */
  inline const bool isHDotPlugged() const
  {
    return hDot->isPlugged();
  }

  /** to set a specified function to compute function hDot
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeHDotFunction(const std::string& , const std::string&);

  /** to compute y = h(q,v,t) using plug-in mechanism
   * \param: double, current time
   */
  void computeH(double);

  /** to compute hDot using plug-in mechanism
   * \param: double, current time
   */
  void computeHDot(double);

  /** to compute the jacobian of h using plug-in mechanism. Index shows which jacobian is computed
   * \param: double, current time
   * \param: unsigned int
   */
  void computeJacH(double, unsigned int);

  /** to compute output
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(double, unsigned int = 0);

  /** to compute p
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  void computeInput(double, unsigned int = 0);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation* : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianRheonomousR* convert(Relation *r);
};

TYPEDEF_SPTR(LagrangianRheonomousR);

#endif
