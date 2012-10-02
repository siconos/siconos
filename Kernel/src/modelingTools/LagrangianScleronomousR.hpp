/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
/*! \file LagrangianScleronomousR.hpp

 */
#ifndef LagrangianScleronomousR_H
#define LagrangianScleronomousR_H

#include "LagrangianR.hpp"

/** Lagrangian (Non Linear) Relations, Scleronomous and Holonomic.
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date Apr 27, 2004
 *
 * Relation with:
 *
 * \f[
 * y = h(q,z) \\
 * \f]
 *
 * \f[
 * \dot y = G0(q,z)\dot q \\
 * \f]
 *
 * \f[
 * p = G0^t(q,z)\lambda
 * \f]
 *
 * with
 *
\f[
G0(q,z) = \nabla_q h(q,z)
\f]
 *
 *  y (or its discrete approximation) is usually stored in y[0]
 *  \dot y (or its discrete approximation) is usually stored in y[1]
 *  higher level can used for storing higher levels of derivatives.
 *
 * G0 and h are connected to plug-in functions.\n
 * The plugin function to compute h(q,z) needs the following parameters:\n
 * --> sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in the interaction\n
 * --> q : pointer to the first element of q \n
 * --> sizeY : size of vector y (ie of the interaction) \n
 * --> [in,out] y : pointer to the first element of y \n
 * --> sizeZ : size of vector z \n
 * --> [in,out] z: pointer to z vector(s) from DS. \n
 * Its signature must be "void plugin(unsigned int, const double*, unsigned int, double*, unsigned int, double*)"\n\n
 * The plugin function to compute G0(q,z), gradient of h according to q, needs the following parameters: \n
 *--> sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in the interaction  \n
 *--> q : pointer to the first element of q  \n
 *--> sizeY : size of vector y (ie of the intercation) \n
 *--> [in,out] G0 : pointer to the first element of G0 (sizeY X sizeDS matrix)\n
 * --> sizeZ : size of vector z \n
 * -->[in,out] z: pointer to z vector(s) from DS.\n
 * Its signature must be "void plugin(unsigned int, const double*, unsigned int, double*, unsigned int, double*)"\n

 *
 */
class LagrangianScleronomousR : public LagrangianR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianScleronomousR);

  /** LagrangianScleronomousR plug-in to compute h(q,z)
  * @param sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in the interaction
  * @param q : pointer to the first element of q
  * @param sizeY : size of vector y (ie of the interaction)
  * @param[in,out] y : pointer to the first element of y
  * @param sizeZ : size of vector z
  * @param[in,out] z: pointer to z vector(s) from DS.
  */
  //  FPtr3 hPtr;

  /** LagrangianScleronomousR plug-in to compute G0(q,z), gradient of h according to q
  * @param sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in the interaction
  * @param q : pointer to the first element of q
  * @param sizeY : size of vector y (ie of the intercation)
  * @param[in,out] G0 : pointer to the first element of G0 (sizeY X sizeDS matrix)
  * @param sizeZ : size of vector z
  * @param[in,out] z: pointer to z vector(s) from DS.
  */
  /** Plugin object for Jacobian of h */
  /** Plugin object for the non linear part of the relative acceleration (derivative of Jacobian of H with
  *respect to the time multiplied by the relative velocity */
  SP::PluggedObject _pluginjqhdot;
  /** Non-linear part of the relative acceleration */
  SP::SiconosVector _NLh2dot;
  /** basic constructor
  \param the sub-type of the relation
  */
  LagrangianScleronomousR(): LagrangianR(RELATION::ScleronomousR)
  {
    ;
  }

  virtual void zeroPlugin();

public:

  /** constructor from xml file
  *  \param relationXML
  */
  LagrangianScleronomousR(SP::RelationXML);

  /** constructor from a set of data
  *  \param string : the name of the plugin to compute h(q,z).\n
  * The signature  of the plugged function must be:
  *  "void pluginH(unsigned int, const double*, unsigned int, double*, unsigned int, double*)"
  *  \param string : the name of the plugin to compute jacobian h according to q.\n
  * The signature  of the plugged function must be:
  *  "void pluginG0(unsigned int, const double*, unsigned int, double*, unsigned int, double*)"
  *  \exception RuntimeException
  */
  LagrangianScleronomousR(const std::string&, const std::string&);

  /** constructor from a set of data used for EventDriven Scheme
  *  \param string : the name of the plugin to compute h(q,z).\n
  * The signature  of the plugged function must be:
  *  "void pluginH(unsigned int, const double*, unsigned int, double*, unsigned int, double*)"
  *  \param string : the name of the plugin to compute jacobian h according to q.\n
  * The signature  of the plugged function must be:
  *  "void pluginG0(unsigned int, const double*, unsigned int, double*, unsigned int, double*)"
  * \param string: the name of the plugin to compute the derivative of H Jacobian with respect to time
  * The signature of the plugged function must be:
  * "void pluginS0(unsigned int, const double*,unsigned int, const double*, unsigned int, double*, unsigned int, double*)"
  *\exception RuntimeException
  */
  LagrangianScleronomousR(const std::string&, const std::string&, const std::string&);

  /** destructor
  */
  virtual ~LagrangianScleronomousR() {};
  /** to get the non-linear part of the relative acceleration */
  inline SP::SiconosVector Nonlinearh2dot()
  {
    return _NLh2dot;
  };

  /** to compute y = h(q,v,t) using plug-in mechanism
  * \param: double, current time
  */
  virtual void computeh(const double time, Interaction& inter);

  /** to compute the jacobian of h using plug-in mechanism. Index shows which jacobian is computed
  * \param: double, current time
  */
  virtual void computeJachq(const double time, Interaction& inter);

  /** to compute the non-linear part of relative accelation using plug-in mechanism
  * \param: double, current time
  * \return: SiconosVector, non linear part
  */
  void computeNonLinearH2dot(const double time, Interaction& inter);

  /** to compute the derivative of H Jacobian with respect to time using plug-in mechanism
  *\param: double, current time
  */
  virtual void computeJachqDot(const double time, Interaction& inter);

  /** to compute output
  *  \param double : current time
  *  \param unsigned int: number of the derivative to compute, optional, default = 0.
  */
  virtual void computeOutput(const double time, Interaction& inter, unsigned int = 0);

  /** to compute p
  *  \param double : current time
  *  \param unsigned int: "derivative" order of lambda used to compute input
  */
  void computeInput(const double time, Interaction& inter, unsigned int = 0);

  const std::string getJachqName() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param Relation* : the relation which must be converted
  * \return a pointer on the relation if it is of the right type, NULL otherwise
  */
  static LagrangianScleronomousR* convert(Relation *r);


  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianScleronomousR)

#endif // LAGRANGIANRELATION_H
