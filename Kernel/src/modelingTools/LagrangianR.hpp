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


/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr2)(unsigned int, double*, unsigned int, double*, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr3)(unsigned int, double*, unsigned int, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr4)(unsigned int, double*, double, unsigned int, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr5bis)(unsigned int, double*, unsigned int, double*, unsigned int, double*, unsigned int, double*);

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
  enum LagrangianRDS  {xfree, z, q0, q1, q2, p0, p1, p2, DSlinkSize};
  // enum LagrangianRVec {vec_xfree, vec_z, vec_q0, vec_q1, vec_q2, vec_p0, vec_p1, vec_p2, vec_workVecSize};
  // enum LagrangianRMat {mat_C, mat_D, mat_F, mat_workMatSize};


protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianR);

  /** Jacobian matrices of \f$y = h(t,q,\dot q,\ldots)\f$ */

  SP::SiconosMatrix _jachlambda;

  /**The Jacobian of the constraints with respect to the generalized coodinates  \f$q\f$
   *  i.e. \f[\nabla^\top_q h(t,q,\dot q,\ldots)\f]
   */
  SP::SimpleMatrix _jachq;

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
   * \param lagType the sub-type of the relation
   */
  LagrangianR(RELATION::SUBTYPES lagType): Relation(RELATION::Lagrangian, lagType) {}

  /** initialize components specific to derived classes.
   * \param inter the interaction using this relation
   * \param DSlink the container of the link to DynamicalSystem attributes
   * \param workV work vectors
   * \param workM work matrices
   */
  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink,
                              VectorOfVectors& workV, VectorOfSMatrices& workM);
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

  /** set Jach[index] to pointer newPtr (pointer link)
  *  \param newPtr the new matrix
  */
  inline void setJachqPtr(SP::SimpleMatrix newPtr)
  {
    _jachq = newPtr ;
  }

  inline SP::SimpleMatrix C() const
  {
    return _jachq;
  }

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction using this relation
   * \param DSlink the container of the link to DynamicalSystem attributes
   * \param workV work vectors
   * \param workM work matrices
  */
  void initialize(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

  /* compute all the H Jacobian 
   * \param time
   * \param inter
   * \param interProp
   */
  virtual void computeJach(double time, Interaction& inter, InteractionProperties& interProp) = 0 ;
  /* compute all the G Jacobian
   * \param time
   * \param inter
   * \param interProp
   */
  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp) = 0 ;

  /** to compute output
   * \param time current time
   * \param inter
   * \param interProp
   *  \param derivativeNumber number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber = 0) = 0;

  /** to compute p
   *  \param time current time
   * \param inter
   * \param interProp
   *  \param level "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0) = 0;

  /** main relation members display
  */
  void display() const;

};
TYPEDEF_SPTR(LagrangianR)
#endif // LAGRANGIANRELATION_H
