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

/*! \file FirstOrderR.h
\brief General interface for relations.
 */

#ifndef FirstOrderR_H
#define FirstOrderR_H

#include "Relation.h"

/** FirstOrder Non Linear Relation.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27, 2004
 *
 *  Relation for First Order Dynamical Systems, with:
 * \f{eqnarray}
 * y &=& h(X,t,\lambda,Z)\\
 * R &=& g(X,t,\lambda,Z)
 * \f}
 *  X, Z, R corresponds to DynamicalSystem variables.
 *  If DS1 and DS2 are involved in the linked Interaction, then X =[x1 x2], Z=[z1 z2] ...
 *
 *  \f$ y \ and \ \lambda \f$ are specific variables of the interaction (see this class for more details).
 *  h and g are plugged on external functions, via plug-in mechanism (see SiconosSharedLibrary).
 *
 * h <=> output
 *
 * g <=> input
 *
 * Operators (and their corresponding plug-in):
     - h: saved in Interaction as y (plug-in: output[0])
     - \f$ \nabla_x h \f$: jacobianH[0] ( output[1] )
     - \f$ \nabla_\lambda h \f$: jacobianH[1] ( output[2] )
     - g: saved in DS as r ( input[0])
     - \f$ \nabla_\lambda g \f$: jacobianG[0] ( input[1] )


     Note: we use a vector for jacobianG while there is only one jacobian. Just for future changes and to allow easy new implementations if some other
     variables are required in g.

 *
 */
template <class T> class FirstOrderR : public Relation
{
public:

  typedef PluggedObject<T, SimpleMatrix> PluggedMatrix;
  typedef boost::shared_ptr<PluggedMatrix> SP_PluggedMatrix;
  typedef T LocalFunc;
  typedef enum DataNames {z, x, r, sizeDataNames};

protected:

  /** Plug-in to compute h(x,t,lambda,z)
   *  @param the size of the vector x.
   *  @param x : the pointer to the first element of the vector x.
   *  @param the size of the vectors y and lambda.
   *  @param[in,out]  a pointer to the first element of the result y
   *  @param the size of the vectors z.
   *  @param[in,out] z : a vector of user-defined parameters.
   */
  T output;

  /** Plug-in to compute g(lambda,t,z)
   *  @param sizeY : the size of the vector y and lambda.
   *  @param lambda : the pointer to the first element of the vector lambda.
   *  @param the size of the vectors R
   *  @param[in,out] : the pointer to the first element of g or its jacobian.
   *  @param the size of the vectors z.
   *  @param[in,out] : a vector of user-defined parameters.
   */
  T input;

  std::vector<SP_PluggedMatrix> JacH;

  std::vector<SP_PluggedMatrix> JacG;

  /** basic constructor
   *  \param the type of the relation
   */
  FirstOrderR(RELATION::SUBTYPES newType): Relation(RELATION::FirstOrder, newType), output(NULL), input(NULL) {}

  /** xml constructor
   *  \param SP::RelationXML : the XML object.
   *  \param the type of the relation
   */
  FirstOrderR(SP::RelationXML relxml, RELATION::SUBTYPES newType): Relation(relxml, RELATION::FirstOrder, newType), output(NULL), input(NULL) {}

  /** To initialize data member: links to DS variables.
   */
  void initDSLinks();

public:

  /** destructor
   */
  virtual ~FirstOrderR() {};

  // -- JacH --

  /** get matrix JacH[index]
   *  \return a SimpleMatrix
   */
  virtual inline const SimpleMatrix getJacH(unsigned int  index = 0) const
  {
    return *(JacH.at(index));
  }

  /** get a pointer on matrix JacH[index]
   *  \return a pointer on a SiconosMatrix
   */
  virtual inline SP::SiconosMatrix getJacHPtr(unsigned int index = 0) const
  {
    return JacH.at(index);
  }

  /** set the value of JacH[index] to newValue (copy)
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in JacH vector
   */
  template <class U> void setJacH(const U& newValue, unsigned int index)
  {
    assert(index < JacH.size() && "FirstOrderR:: setJacH(mat,index), index out of range. Maybe you do not set the sub-type of the relation?");
    // resize the local matrix, else exception. Note FP: Set resize in SimpleMatrix::operator= ??
    if (JacH[index])
      JacH[index]->resize(newValue.size(0), newValue.size(1));
    setObject<PluggedMatrix, SP_PluggedMatrix, U>(JacH[index], newValue);
  }

  /** set JacH[index] to pointer newPtr (pointer link)
   *  \param SP::SiconosMatrix  newPtr
   *  \param unsigned int: index position in JacH vector
   */
  inline void setJacHPtr(SP_PluggedMatrix newPtr, unsigned int index = 0)
  {
    JacH.at(index) = newPtr ;
  }

  // -- JacG --

  /** get matrix JacG[index]
   *  \return a SimpleMatrix
   */
  inline const SimpleMatrix getJacG(unsigned int  index = 0) const
  {
    return *(JacG.at(index));
  }

  /** get a pointer on matrix JacG[index]
   *  \return a pointer on a SiconosMatrix
   */
  virtual inline SP::SiconosMatrix getJacGPtr(unsigned int index = 0) const
  {
    return JacG.at(index);
  }

  /** set the value of JacG[index] to newValue (copy)
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in JacG vector
   */
  template <class U> void setJacG(const U& newValue, unsigned int index)
  {
    assert(index < JacG.size() && "FirstOrderR:: setJacG(mat,index), index out of range. Maybe you do not set the sub-type of the relation?");
    if (JacG[index]) JacG[index]->resize(newValue.size(0), newValue.size(1));
    setObject<PluggedMatrix, SP_PluggedMatrix, U>(JacG[index], newValue);
  };

  /** set JacG[index] to pointer newPtr (pointer link)
   *  \param SP::SiconosMatrix  newPtr
   *  \param unsigned int: index position in JacG vector
   */
  inline void setJacGPtr(SP_PluggedMatrix newPtr, unsigned int index = 0)
  {
    JacG.at(index) = newPtr ;
  }

  /** To get the name of JacH[i] plugin
   *  \return a string
   */
  const std::string getJacHName(unsigned int i) const
  {
    return JacH[i]->getPluginName();
  }

  /** To get the name of JacG[i] plugin
   *  \return a string
   */
  const std::string getJacGName(unsigned int i) const
  {
    return JacG[i]->getPluginName();
  }

  /** true if JacH[i] is plugged
   *  \return a bool
   */
  const bool isJacHPlugged(unsigned int i) const
  {
    return JacH[i]->isPlugged();
  }

  /** true if JacG[i] is plugged
   *  \return a bool
   */
  const bool isJacGPlugged(unsigned int i) const
  {
    return JacG[i]->isPlugged();
  }

  /** Gets the number of computed jacobians for h
      \return an unsigned int.
  */
  inline unsigned int getNumberOfJacobiansForH() const
  {
    return JacH.size();
  }

  /** Gets the number of computed jacobian for g
      \return an unsigned int.
  */
  inline unsigned int getNumberOfJacobiansForG() const
  {
    return JacG.size();
  }

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
  void setComputeJacobianHFunction(const std::string&, const std::string&, unsigned int = 0);

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

  /** initialize the relation (check sizes, memory allocation ...)
      \param SP to Interaction: the interaction that owns this relation
   */
  void initialize(SP::Interaction);

  /** default function to compute h
   *  \param double : current time
   */
  void computeH(double) = 0;

  /** default function to compute g
   *  \param double : current time
   */
  void computeG(double) = 0;

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  void computeJacH(double, unsigned int);

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  void computeJacG(double, unsigned int);

  /** default function to compute y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(double, unsigned int = 0) = 0;

  /** default function to compute r
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  void computeInput(double, unsigned int = 0) = 0;

  /** main relation members display
   */
  void display() const;
};

#endif
