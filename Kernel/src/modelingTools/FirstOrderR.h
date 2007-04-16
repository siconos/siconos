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

/*! \file FirstOrderR.h
\brief General interface for relations.
 */

#ifndef FirstOrderR_H
#define FirstOrderR_H

#include "Relation.h"

/** Pointer to function for plug-in for operators related to output and its gradients.*/
typedef void (*OutPtr)(unsigned int, const double*, double, unsigned int, const double*, double*, unsigned int, double*);

/** Pointer to function for plug-in for operators related to input and its gradients.*/
typedef void (*InPtr)(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*);

/** FirstOrder Non Linear Relation.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) Apr 27, 2004
 *
 *  Relation for First Order Dynamical Systems, with:
 * \f{eqnarray}
 * y &=& h(X,t,\lambda,Z)\\
 * R &=& g(\lambda,t,Z)
 * \f}
 *  X, Z, R corresponds to DynamicalSystem variables.
 *  If DS1 and DS2 are involved in the linked Interaction, then X =[x1 x2], Z=[z1 z2] ...
 *  \f$ y \ and \ \lambda \f$ are specific variables of the interaction (see this class for more details).
 *  h and g are plugged on external functions, via plug-in mechanism (see SiconosSharedLibrary).
 *
 * h <=> output
 *
 * g <=> input
 *
 * Operator (and its corresponding plug-in):
     - h: saved in Interaction as y (plug-in: output[0])
     - \f$ \nabla_x h \f$: jacobianH[0] ( output[1] )
     - \f$ \nabla_\lambda h \f$: jacobianH[1] ( output[2] )
     - g: saved in DS as R ( input[0])
     - \f$ \nabla_\lambda g \f$: jacobianG[0] ( input[1] )


     Note: we use a vector for jacobianG while there is only one jacobian. Just for future changes and to allow easy new implementations if some other
     variables are required in g.

 *
 */
class FirstOrderR : public Relation
{

protected:

  /** type of the FirstOrderR */
  std::string  firstOrderType;

  /** matrices of gradients of h. (jacobianH[0]: gradient according to x, jacobianH[1]: gradient according to lambda) */
  VectorOfMatrices jacobianH;

  /** matrices of gradients of g. (jacobianG[0]: gradient according to x, jacobianG[1]: gradient according to lambda) */
  VectorOfMatrices jacobianG;

  /** Plug-in to compute h(x,t,lambda,z) (index = 0), \f$ \nabla_x h(x,t,lambda,z)\f$ (index = 1) or \f$ \nabla_\lambda h(x,t,lambda,z)\f$ (index = 2)
   *  @param the size of the vector x.
   *  @param x : the pointer to the first element of the vector x.
   *  @param time : current time.
   *  @param the size of the vectors y and lambda.
   *  @param lambda : the pointer to the first element of the vector lambda.
   *  @param[in,out]  a pointer to the first element of the result (y or its jacobian)
   *  @param the size of the vectors z.
   *  @param[in,out] z : a vector of user-defined parameters.
   */
  std::vector<OutPtr> output;

  /** Plug-in to compute g(lambda,t,z) (index = 0), \f$ \nabla_X g(lambda,t,z)\f$ (index = 1) or \f$ \nabla_\lambda g(lambda,t,z)\f$ (index = 2).
   *  @param sizeY : the size of the vector y and lambda.
   *  @param lambda : the pointer to the first element of the vector lambda.
   *  @param time : current time.
   *  @param the size of the vectors R
   *  @param[in,out] : the pointer to the first element of g or its jacobian.
   *  @param the size of the vectors z.
   *  @param[in,out] : a vector of user-defined parameters.
   */
  std::vector<InPtr> input;

  /** protected function used to initilized isPlugged map
      \param : a bool, value for all flags.
  */
  virtual void initPluginFlags(bool);

  /** protected function used to initilized isAllocated map
      \param : a bool, value for all flags.
  */
  virtual void initAllocationFlags(bool);

  /** default constructor
   *  \param a string that gives the type of the relation (optional)
   */
  FirstOrderR(const std::string& = "R");

public:

  /** xml constructor
   *  \param FirstOrderRXML* : the XML object.
   *  \param a string that gives the type of the relation (optional)
   */
  FirstOrderR(RelationXML*, const std::string& = "R");

  /** data constructor
   *  \param a string with computeOutput function name.
   *  \param a string with computeInput function name.
   */
  FirstOrderR(const std::string&, const std::string&);

  /** destructor
   */
  virtual ~FirstOrderR();

  /** To get the type of the FirstOrderR
   *  \return string : the type of the FirstOrderR
   */
  inline const std::string  getFirstOrderRelationType() const
  {
    return firstOrderType;
  }

  /** initialize the relation (check sizes, memory allocation ...)
   */
  virtual void initialize();

  // -- jacobianH --

  /** get the vector of matrices jacobianH
   *  \return vector<SiconosMatrix*>
   */
  inline VectorOfMatrices getJacobianHVector() const
  {
    return jacobianH;
  }

  /** set the values of jacobianH[i]: no copy but links with pointers!!
   *  \param a VectorOfMatrices.
   */
  void setJacobianHVector(const VectorOfMatrices&);

  /** get matrix jacobianH[index]
   *  \return a SimpleMatrix
   */
  inline const SimpleMatrix getJacobianH(unsigned int  index = 0) const
  {
    return *(jacobianH[index]);
  }

  /** get a pointer on matrix jacobianH[index]
   *  \return a pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianHPtr(unsigned int index = 0) const
  {
    return jacobianH[index];
  }

  /** set the value of jacobianH[index] to newValue (copy)
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in jacobianH vector
   */
  void setJacobianH(const SiconosMatrix&, unsigned int = 0);

  /** set jacobianH[index] to pointer newPtr (pointer link)
   *  \param SiconosMatrix * newPtr
   *  \param unsigned int: index position in jacobianH vector
   */
  void setJacobianHPtr(SiconosMatrix *newPtr, unsigned int = 0);

  // -- jacobianG --

  /** get the vector of matrices jacobianG
   *  \return vector<SiconosMatrix*>
   */
  inline VectorOfMatrices getJacobianGVector() const
  {
    return jacobianG;
  }

  /** set the values of jacobianG[i]: no copy but links with pointers!!
   *  \param a VectorOfMatrices.
   */
  void setJacobianGVector(const VectorOfMatrices&);

  /** get matrix jacobianG[index]
   *  \return a SimpleMatrix
   */
  inline const SimpleMatrix getJacobianG(unsigned int  index = 0) const
  {
    return *(jacobianG[index]);
  }

  /** get a pointer on matrix jacobianG[index]
   *  \return a pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianGPtr(unsigned int index = 0) const
  {
    return jacobianG[index];
  }

  /** set the value of jacobianG[index] to newValue (copy)
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in jacobianG vector
   */
  void setJacobianG(const SiconosMatrix&, unsigned int = 0);

  /** set jacobianG[index] to pointer newPtr (pointer link)
   *  \param SiconosMatrix * newPtr
   *  \param unsigned int: index position in jacobianG vector
   */
  void setJacobianGPtr(SiconosMatrix *newPtr, unsigned int = 0);

  /** To set a plug-in function to compute output function h
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeHFunction(const std::string&, const std::string&);

  /** To set a plug-in function to compute jacobianH
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  void setComputeJacobianHFunction(const std::string&, const std::string&, unsigned int);

  /** To set a plug-in function to compute input function g
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeGFunction(const std::string&, const std::string&);

  /** To set a plug-in function to compute the jacobian according to x of the input
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  void setComputeJacobianGFunction(const std::string&, const std::string&, unsigned int = 0);

  /** default function to compute y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double, unsigned int = 0);

  /** default function to compute y for the free state
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeFreeOutput(double, unsigned int = 0);

  /** default function to compute r
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double, unsigned int);

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJacobianH(double, unsigned int);

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  virtual void computeJacobianG(double, unsigned int = 0);

  /** main relation members display
   */
  virtual void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static FirstOrderR* convert(Relation *r);
};

#endif
