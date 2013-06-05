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
/*! \file LagrangianR.hpp

 */
#ifndef LAGRANGIANRELATION_H
#define LAGRANGIANRELATION_H

#include "Relation.hpp"
#include "Interaction.hpp"

class DynamicalSystem;
class RelationXML;
class SimpleMatrix;
class SiconosVector;

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr2)(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr3)(unsigned int, const double*, unsigned int, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr4)(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr5bis)(unsigned int, const double*, unsigned int, const double*, unsigned int, double*, unsigned int, double*);

/** Lagrangian (Non Linear) Relation (generic interface)
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date Apr 27, 2004
 *
 * Relations for Lagrangian Dynamical Systems.
 * This class is only an interface for specific (Linear, scleronomic, rheomomic  ...)
 * Lagrangian Relations (see derived classes).
 *
 *  Class name = type+subType.
 *
 * If \f$y = h(t,q,\dot q,\ldots)\f$ describes the constraint (the relation) , all the gradients of h
 * are handled by the following SiconosMatrix and SiconosVector objects.
 *
 * <ul>
 * <li> The Jacobian of the constraints with respect to the coodinates  \f$q\f$
 * i.e. \f[\nabla^T_q h(t,q,\dot q,\ldots)\f]  is stored in  SP::SiconosMatrix _jachq .
 *
 * This Jacobian is mainly used for Newton linearization and to compute the time-derivative of the constraint \f$y = h(q,\ldots)\f$ that is
 *  \f[\dot y (t) = \nabla^T_q h(t,q,\dot q,\ldots) (q) \dot q +\ldots\f]
 * This object can also store
 * more general linearized part of the gap function. If \f$y=h(q)\f$ models a gap function, then the time--derivative
 * can be generically  written as
 * \f[\dot y (t) = H(q,\ldots) \dot q  +\ldots. \f]
 * The matrix \f$H(q,\ldots) \f$ is also stored in   SP::SiconosMatrix _jachq </li>
 *
 * <li> The Jacobian of the constraints with respect to the generalized velocities  \f$\dot q\f$
 *  i.e. \f[\nabla^\top_{\dot q} h(t,q,\dot q,\ldots)\f] is stored in  SP::SiconosMatrix _jachqDot </li>
 *
 * <li>The time-derivative of Jacobian of the constraints with respect to the generalized coordinates  \f$ q\f$
 *  i.e. \f[\frac{d}{dt} \nabla^\top_{q} h(t,q,\dot q,\ldots).\f]. This value is useful to compute the second-order
 * time--derivative of the constraints with respect to time.</li>
 *
 * </ul>
 *
 * In corresponding derived classes, h and Jacobians are connected to plug-in functions (user-defined).
 *
 */
class LagrangianR : public Relation
{
public:

  enum DataNames {free, z, q0, q1, q2, p0, p1, p2, sizeDataNames};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianR);

  /** Jacobian matrices of \f$y = h(t,q,\dot q,\ldots)\f$ */

  /**The Jacobian of the constraints with respect to the generalized coodinates  \f$q\f$
   *  i.e. \f[\nabla^\top_q h(t,q,\dot q,\ldots)\f]
   */
  SP::SiconosMatrix _jachq;

  /**The Jacobian of the constraints with respect to the generalized velocities  \f$\dot q\f$
   *  i.e. \f[\nabla^\top_{\dot q} h(t,q,\dot q,\ldots)\f]
   */
  SP::SiconosMatrix _jachqDot;

  /**The time-derivative of Jacobian of the constraints with respect
     to the generalized coordinates  \f$ q\f$
   * i.e. \f[\frac{d}{dt} \nabla^\top_{ q} h(t,q,\dot q,\ldots).\f]
   * This value is useful to compute the second-order
   * time--derivative of the constraints with respect to time.
   */
  SP::SiconosMatrix _dotjachq;

  SP::PluggedObject _pluginJachq;

  /** basic constructor
  \param the sub-type of the relation
  */
  LagrangianR(RELATION::SUBTYPES lagType): Relation(RELATION::Lagrangian, lagType) {}

  /** constructor from xml file
  *  \param relationXML
  *  \param std::string: relation subType
  */
  LagrangianR(SP::RelationXML relxml, RELATION::SUBTYPES newSubType): Relation(relxml, RELATION::Lagrangian, newSubType) {}

  /** initialize components specific to derived classes.
  */
  virtual void initComponents(Interaction& inter);
  virtual void zeroPlugin();

public:

  /** destructor
  */
  virtual ~LagrangianR() {};

  // -- Jach --

  /** get matrix Jach[index]
  *  \return a SimpleMatrix
  inline const SimpleMatrix getJach(unsigned int  index = 0) const { return *(Jach.at(index)); }
  */

  /** get a pointer on matrix Jach[index]
  *  \return a pointer on a SiconosMatrix
  */
  inline SP::SiconosMatrix jachq() const
  {
    return _jachq;
  }
  inline SP::SiconosMatrix jachqDot() const
  {
    return _jachqDot;
  }
  inline SP::SiconosMatrix dotJachq() const
  {
    return _dotjachq;
  }
  inline SP::SiconosMatrix jachlambda() const
  {
    return _jachlambda;
  }

  /** set the value of Jach[index] to newValue (copy)
  *  \param SiconosMatrix newValue
  *  \param unsigned int: index position in Jach vector

  template <class U> void setJach(const U& newValue, unsigned int index = 0)
  {
  assert(index>=Jach.size()&&"LagrangianR:: setJach(mat,index), index out of range. Maybe you do not set the sub-type of the relation?");

  if(Jach[index]) Jach[index]->resize(newValue.size(0), newValue.size(1));
  setObject<PluggedMatrix,SP_PluggedMatrix,U>(Jach[index],newValue);
  }
  */
  /** set Jach[index] to pointer newPtr (pointer link)
  *  \param SP::SiconosMatrix  newPtr
  *  \param unsigned int: index position in Jach vector
  */
  inline void setJachqPtr(SP::SiconosMatrix newPtr)
  {
    _jachq = newPtr ;
  }

  /** To get the name of Jach[i] plugin
  *  \return a std::string
  const std::string getJachName(unsigned int i) const {return Jach[i]->getPluginName();}
  */

  /** true if Jach[i] is plugged
  *  \return a bool
  const bool isJachPlugged(unsigned int i) const {return Jach[i]->isPlugged();}
  */

  /** Gets the number of computed jacobians for h
  \return an unsigned int.
  inline unsigned int numberOfJacobiansForH() const { return Jach.size();}
  */

  inline SP::SiconosMatrix C() const
  {
    return _jachq;
  }

  /** initialize the relation (check sizes, memory allocation ...)
  \param SP to Interaction: the interaction that owns this relation
  */
  void initialize(Interaction& inter);

  /** to compute y = h(t,q,v,z) using plug-in mechanism
   * should be used as less as possible to avoid side--effects
   * prefer computeh(const double time, Interaction& inter,
                     SP::BlockVector q, SP::BlockVector v, SP::BlockVector z)
   * \param time  current time
   * \param inter interaction that owns the relation
   */
  virtual void computeh(const double time, Interaction& inter);

  /** to compute y = h(t,q,v,z) using plug-in mechanism
  * \param time current time
  * \param inter interaction that owns the relation
  * \param q the BlockVector of coordinates
  * \param v the BlockVector of velocities
  * \param z the BlockVector of parameters
  */
  void computeh(const double time, Interaction& inter,
                SP::BlockVector q, SP::BlockVector v, SP::BlockVector z);

  // void computeh(const double time, Interaction& inter,
  //               SP::BlockVector q, SP::BlockVector v,
  //               SP::BLockVector lambda, SP::BlockVector z
  //               SP::SiconosVector y);
  /** default function to compute jacobianH
  *  \param double : current time
  *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)

  void computeJachx(double);*/
  virtual void computeJachlambda(const double time, Interaction& inter)
  {
    ;
  }
  virtual void computeJachq(const double time, Interaction& inter)
  {
    ;
  }
  virtual void computeJachqDot(const double time, Interaction& inter)
  {
    ;
  }
  virtual void computeDotJachq(const double time, Interaction& inter)
  {
    ;
  }

  /** to compute hDot using plug-in mechanism
   * using plug-in mechanism with the data vector of the interaction
   * should be used as less as possible to avoid side--effects
   * prefer computehDot(const double time, Interaction& inter, SP::BlockVector q, SP::BlockVector z)
   * \param time  current time
   * \param inter interaction that owns the relation
   */
  virtual void computehDot(const double time, Interaction& inter)
  {
    ;
  }

  void computeJacglambda(const double time, Interaction& inter)
  {
    ;
  }
  void computeJacgq(const double time, Interaction& inter)
  {
    ;
  }
  void computeJacgqDot(const double time, Interaction& inter)
  {
    ;
  }
  /* compute all the H Jacobian */
  virtual void computeJach(const double time, Interaction& inter)
  {
    computeJachq(time, inter);
    computeJachqDot(time, inter);
    computeDotJachq(time, inter);
    computeJachlambda(time, inter);
    computehDot(time,inter);
  }
  /* compute all the G Jacobian */
  virtual void computeJacg(const double time, Interaction& inter)
  {
    computeJacgq(time, inter);
    computeJacgqDot(time, inter);
    computeJacglambda(time, inter);
  }


  /** to compute output
  *  \param double : current time
  *  \param unsigned int: number of the derivative to compute, optional, default = 0.
  */
  virtual void computeOutput(const double time, Interaction& inter, unsigned int = 0) = 0;

  /** to compute p
  *  \param double : current time
  *  \param unsigned int: "derivative" order of lambda used to compute input
  */
  virtual void computeInput(const double time, Interaction& inter, unsigned int = 0) = 0;

  /** copy the data of the Relation to the XML tree
  */
  void saveRelationToXML() const;

  /** main relation members display
  */
  void display() const;

};
TYPEDEF_SPTR(LagrangianR)
#endif // LAGRANGIANRELATION_H
