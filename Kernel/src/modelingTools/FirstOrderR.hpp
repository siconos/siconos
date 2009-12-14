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

#include "Relation.hpp"

/** Pointer to function for plug-in for operators related to output and its gradients.*/
typedef void (*OutPtr)(unsigned int, const double*, double, unsigned int, const double*, double*, unsigned int, double*);

/** Pointer to function for plug-in for operators related to input and its gradients.*/
typedef void (*InPtr)(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*);

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
class FirstOrderR : public Relation
{
public:

  enum DataNames {z, x, r, g_alpha, ds_xp, sizeDataNames};

protected:



  SP::SiconosMatrix JacXH;

  SP::SiconosMatrix JacLG;

  /** basic constructor
   *  \param the type of the relation
   */
  FirstOrderR(RELATION::SUBTYPES newType): Relation(RELATION::FirstOrder, newType) {}

  /** xml constructor
   *  \param SP::RelationXML : the XML object.
   *  \param the type of the relation
   */
  FirstOrderR(SP::RelationXML relxml, RELATION::SUBTYPES newType): Relation(relxml, RELATION::FirstOrder, newType) {}

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

  virtual inline const SimpleMatrix getJacXH() const { return *(JacH.at(index)); }
  */
  /** get a pointer on matrix JacH[index]
   *  \return a pointer on a SiconosMatrix
   */
  virtual  SP::SiconosMatrix jacXH() const
  {
    return JacXH;
  }

  /** set the value of JacH[index] to newValue (copy)
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in JacH vector

  void setJacobianH(const SiconosMatrix&, unsigned int = 0);
  */

  /** set JacH[index] to pointer newPtr (pointer link)
   *  \param SP::SiconosMatrix  newPtr
   *  \param unsigned int: index position in JacH vector
   */
  inline void setJacXHPtr(SP::SiconosMatrix newPtr)
  {
    JacXH = newPtr ;
  }
  inline void setJacLHPtr(SP::SiconosMatrix newPtr)
  {
    JacLH = newPtr ;
  }

  // -- JacG --

  /** get matrix JacG[index]
   *  \return a SimpleMatrix

  inline const SimpleMatrix getJacG(unsigned int  index = 0) const { return *(JacG.at(index)); }
  */
  /** get a pointer on matrix JacG[index]
   *  \return a pointer on a SiconosMatrix
   */
  virtual SP::SiconosMatrix jacLG() const
  {
    return JacLG;
  }

  /** set the value of JacG[index] to newValue (copy)
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in JacG vector

  void setJacG(const U& newValue, unsigned int index )
    {
      assert(index<JacG.size()&&"FirstOrderR:: setJacG(mat,index), index out of range. Maybe you do not set the sub-type of the relation?");
      if(JacG[index]) JacG[index]->resize(newValue.size(0), newValue.size(1));
      setObject<PluggedMatrix,SP_PluggedMatrix,U>(JacG[index],newValue);
    };
  */

  /** set JacG[index] to pointer newPtr (pointer link)
   *  \param SP::SiconosMatrix  newPtr
   *  \param unsigned int: index position in JacG vector
   */
  inline void setJacLGPtr(SP::SiconosMatrix newPtr)
  {
    JacLG = newPtr ;
  }

  /** To get the name of JacH[i] plugin
   *  \return a string
  const std::string getJacHName(unsigned int i) const {return JacH[i]->getPluginName();}
   */

  /** To get the name of JacG[i] plugin
   *  \return a string
  const std::string getJacGName(unsigned int i) const {return JacG[i]->getPluginName();}
   */

  /** true if JacH[i] is plugged
   *  \return a bool
   */

  /** initialize the relation (check sizes, memory allocation ...)
      \param SP to Interaction: the interaction that owns this relation
   */
  virtual void initialize(SP::Interaction);

  /** default function to compute h
   *  \param double : current time
   */
  void computeh(double) = 0;

  /** default function to compute g
   *  \param double : current time
   */
  void computeg(double) = 0;

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJacXH(double);
  virtual void computeJacLH(double);
  virtual void computeJacH(double t)
  {
    computeJacXH(t);
    computeJacLH(t);
  }


  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  void computeJacLG(double);

  virtual void computeJacG(double t)
  {
    computeJacLG(t);
  }


  /** main relation members display
   */

  /**
  * return a SP on the C matrix.
  * The matrix C in the linear case, else it returns Jacobian of the output with respect to x.
  *
  */
  SP::SiconosMatrix getCPtr()
  {
    return jacXH();
  }
  /**
   * return a SP on the D matrix.
   * The matrix D in the linear case, else it returns Jacobian of the output with respect to lambda.
   */
  SP::SiconosMatrix getDPtr()
  {
    return jacLH();
  }
  /**
   * return a SP on the B matrix.
   * The matrix B in the linear case, else it returns Jacobian of the input with respect to lambda.
   */
  SP::SiconosMatrix getBPtr()
  {
    return jacLG();
  }


  void display() const;
};
TYPEDEF_SPTR(FirstOrderR);

#endif
