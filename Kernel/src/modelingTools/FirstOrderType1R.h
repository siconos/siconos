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

/*! \file FirstOrderType1.h
\brief non linear relations, with y depending on dynamical systems state and r on lambda.
 */

#ifndef FirstOrderType1R_H
#define FirstOrderType1R_H

#include "FirstOrderR.h"

/** Pointer to function for plug-in for operators related to output and its gradients.*/
typedef void (*Type1Ptr)(unsigned int, const double*, unsigned int, double*, unsigned int, double*);

/** FirstOrder Non Linear Relation.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27, 2004
 *
 * Derived from FirstOrderR - See this class for more comments.
 *
 *  Relation for First Order Dynamical Systems, with:
 * \f{eqnarray}
 * y &=& h(X,Z)\\
 * r &=& g(\lambda,Z)
 * \f}
 *
 * Operators (and their corresponding plug-in):
     - h: saved in Interaction as y (plug-in: output[0])
     - \f$ \nabla_x h \f$: jacobianH[0] ( output[1] )
     - g: saved in DS as r ( input[0])
     - \f$ \nabla_\lambda g \f$: jacobianG[0] ( input[1] )
 *
 */
class FirstOrderType1R : public FirstOrderR
{

protected:

  /** Plug-in to compute h(x,t,lambda,z)
   *  @param the size of the vector x.
   *  @param x : the pointer to the first element of the vector x.
   *  @param the size of the vectors y and lambda.
   *  @param[in,out]  a pointer to the first element of the result y
   *  @param the size of the vectors z.
   *  @param[in,out] z : a vector of user-defined parameters.
   */
  Type1Ptr output;

  /** Plug-in to compute \f$ \nabla_x h(x,t,lambda,z)\f$
   *  @param the size of the vector x.
   *  @param x : the pointer to the first element of the vector x.
   *  @param the size of the vectors y and lambda.
   *  @param[in,out]  a pointer to the first element of the result, the jacobian
   *  @param the size of the vectors z.
   *  @param[in,out] z : a vector of user-defined parameters.
   */
  Type1Ptr jXOutput;

  /** Plug-in to compute g(lambda,t,z)
   *  @param sizeY : the size of the vector y and lambda.
   *  @param lambda : the pointer to the first element of the vector lambda.
   *  @param the size of the vectors R
   *  @param[in,out] : the pointer to the first element of g or its jacobian.
   *  @param the size of the vectors z.
   *  @param[in,out] : a vector of user-defined parameters.
   */
  Type1Ptr input;

  /** Plug-in to compute \f$ \nabla_\lambda g(lambda,t,z)\f$
   *  @param sizeY : the size of the vector y and lambda.
   *  @param lambda : the pointer to the first element of the vector lambda.
   *  @param the size of the vectors R
   *  @param[in,out] : the pointer to the first element of the jacobian.
   *  @param the size of the vectors z.
   *  @param[in,out] : a vector of user-defined parameters.
   */
  Type1Ptr jLInput;

  /** protected function used to initilized isPlugged map
      \param : a bool, value for all flags.
  */
  void initPluginFlags(bool);

private:

  /** default constructor
   */
  FirstOrderType1R();

public:

  /** xml constructor
   *  \param FirstOrderType1RXML* : the XML object.
   */
  FirstOrderType1R(RelationXMLSPtr);

  /** data constructor
   *  \param a string with computeOutput function name.
   *  \param a string with computeInput function name.
   */
  FirstOrderType1R(const std::string&, const std::string&);

  /** data constructor
   *  \param a string with computeOutput function name.
   *  \param a string with computeInput function name.
   *  \param a string: name of the function to compute the jacobian of h according to x
   *  \param a string: name of the function to compute the jacobian of g according to lambda
   */
  FirstOrderType1R(const std::string&, const std::string&, const std::string&, const std::string&);

  /** destructor
   */
  ~FirstOrderType1R();

  /** initialize the relation (check sizes, memory allocation ...)
   */
  void initialize();

  /** To set a plug-in function to compute output function h
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeHFunction(const std::string&, const std::string&);

  /** To set a plug-in function to compute jacobian of h according to x
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \param index: not used, set any value.
   */
  void setComputeJacobianHFunction(const std::string&, const std::string&, unsigned int = 0);

  /** To set a plug-in function to compute input function g
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeGFunction(const std::string&, const std::string&);

  /** To set a plug-in function to compute the jacobian of g according to lambda
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \param index not used.
   */
  void setComputeJacobianGFunction(const std::string&, const std::string&, unsigned int = 0);

  /** default function to compute y
   *  \param double: not used
   *  \param unsigned int: not used
   */
  void computeOutput(double, unsigned int = 0);

  /** default function to compute r
   *  \param double : not used
   *  \param unsigned int: not used
   */
  void computeInput(double, unsigned int = 0);

  /** default function to compute jacobianH
   *  \param double : not used
   *  \param not used
   */
  void computeJacobianH(double, unsigned int);

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  void computeJacobianG(double, unsigned int);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static FirstOrderType1R* convert(Relation *r);
};

#endif
