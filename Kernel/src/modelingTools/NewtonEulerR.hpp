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
/*! \file NewtonEulerR.hpp

 */
#ifndef NEWTONEULERRELATION_H
#define NEWTONEULERRELATION_H

#include "Relation.hpp"
#include "Interaction.hpp"

class DynamicalSystem;
class RelationXML;
class SimpleMatrix;
class SiconosVector;

/** NewtonEuler (Non Linear) Relation (generic interface)
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date Apr 27, 2004
 *
 * Relations for NewtonEuler Dynamical Systems. This class is only an interface for specific (Linear, Scleronomous ...)
 * NewtonEuler Relations (see derived classes).
 *
 *  Class name = type+subType.
 *
 * If y = h(...), all the gradients of are handled by G object.
 * For example, G[0] = \f$ \nabla_q h(q,...) \f$.
 *
 * In corresponding derived classes, h and Gi are connected to plug-in functions (user-defined).
 * For more details, see the DevNotes.pdf, chapter NewtonEuler.
 */

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr2)(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr3)(unsigned int, const double*, unsigned int, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr4)(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*);


class NewtonEulerR : public Relation
{
public:

  enum DataNames {free, z, q0, q1, q2, p0, p1, p2, velo, deltaq, sizeDataNames};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEulerR);

  /** Jacobian matrices of H */
  SP::SimpleMatrix _jachq;
  //proj_with_q SP::SimpleMatrix _jachqProj;
  SP::SiconosMatrix _jachqDot;
  SP::SiconosMatrix _jacglambda;

  /**vector e*/
  SP::SiconosVector _e;
  /*Used for the projection formulation*/

  /**vector of contact forces, ie: _contactForce = B \lambda. Useful for the end user.*/
  SP::SiconosVector _contactForce;

  /**updated in computeJachqT:
  In the case of the bilateral constrains, it is _jachq._T.
  In the case of a local frame, _jachqT is built from the geometrical datas(local frame, point of contact).*/
  SP::SiconosMatrix _jachqT;

  /** basic constructor
  \param the sub-type of the relation
  */
  NewtonEulerR(RELATION::SUBTYPES lagType): Relation(RELATION::NewtonEuler, lagType) {}

  /** constructor from xml file
  *  \param relationXML
  *  \param string: relation subType
  */
  NewtonEulerR(SP::RelationXML relxml, RELATION::SUBTYPES newSubType): Relation(relxml, RELATION::NewtonEuler, newSubType) {}

  /** initialize components specific to derived classes.
  */
  virtual void initComponents(Interaction& inter);

public:
  NewtonEulerR(): Relation(RELATION::NewtonEuler, RELATION::NonLinearR) {}

  /** destructor
  */
  virtual ~NewtonEulerR() {};

  // -- Jach --

  /** get a pointer on matrix Jach[index]
  *  \return a pointer on a SiconosMatrix
  */
  inline SP::SimpleMatrix jachq() const
  {
    return _jachq;
  }
  //proj_with_q  inline SP::SimpleMatrix jachqProj() const { return _jachqProj; }
  inline void setJachq(SP::SimpleMatrix newJachq)
  {
    _jachq = newJachq;
  }
  inline SP::SiconosMatrix jacQDotH() const
  {
    return _jachqDot;
  }
  inline SP::SiconosMatrix jachlambda() const
  {
    return _jachlambda;
  }
  inline SP::SiconosMatrix jacglambda() const
  {
    return _jacglambda;
  }
  inline void setE(SP::SiconosVector newE)
  {
    _e = newE;
  }

  inline SP::SiconosMatrix jachqT() const
  {
    return _jachqT;
  }
  inline void setJachqT(SP::SiconosMatrix newJachqT)
  {
    _jachqT = newJachqT;
  }

  /** set Jach[index] to pointer newPtr (pointer link)
  *  \param SP::SiconosMatrix  newPtr
  *  \param unsigned int: index position in Jach vector
  */
  inline void setJachqPtr(SP::SimpleMatrix newPtr)
  {
    _jachq = newPtr ;
  }

  /** initialize the relation (check sizes, memory allocation ...)
  \param SP to Interaction: the interaction that owns this relation
  */
  void initialize(Interaction& inter);

  /** to compute y = h(q,v,t) using plug-in mechanism
  * \param: double, current time
  */
  virtual void computeh(const double time, Interaction& inter);

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
  virtual void computeJacglambda(const double time, Interaction& inter)
  {
    ;
  }
  virtual void computeJacgq(const double time, Interaction& inter)
  {
    ;
  }
  virtual void computeJacgqDot(const double time, Interaction& inter)
  {
    ;
  }
  /*default implementation consists in multiplying jachq and T*/
  virtual void computeJachqT(Interaction& inter);
  /* compute all the H Jacobian */
  virtual void computeJach(const double time, Interaction& inter)
  {
    computeJachq(time, inter);
    computeJachqDot(time, inter);
    computeJachlambda(time, inter);
    computeJachqT(inter);
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
  virtual void computeOutput(const double time, Interaction& inter, unsigned int = 0) ;

  /** to compute p
  *  \param double : current time
  *  \param unsigned int: "derivative" order of lambda used to compute input
  */
  virtual void computeInput(const double time, Interaction& inter, unsigned int = 0) ;

  /** copy the data of the Relation to the XML tree
  */
  void saveRelationToXML() const;

  /**
  * return a SP on the C matrix.
  * The matrix C in the linear case, else it returns Jacobian of the output with respect to x.
  *
  */
  virtual inline SP::SiconosMatrix C() const
  {
    return _jachq;
  }
  /**
  * return a SP on the D matrix.
  * The matrix D in the linear case, else it returns Jacobian of the output with respect to lambda.
  */
  virtual inline SP::SiconosMatrix D() const
  {
    return _jachlambda;
  }
  /**
  * return a SP on the B matrix.
  * The matrix B in the linear case, else it returns Jacobian of the input with respect to lambda.
  */
  virtual inline SP::SiconosMatrix B() const
  {
    return _jacglambda;
  }
  /** A buffer containing the forces due to this.
  It is an output unused for the computation.
  */
  inline SP::SiconosVector contactForce() const
  {
    return _contactForce;
  };

  /** return a block vector containing ths dynamical system's dof.*/
  //  SP::BlockVector getq(){return data[q0];}

  //  SP::BlockVector q(){return data[q0];}

  //  SP::BlockVector getDeltaq(){return data[deltaq];}

  ACCEPT_STD_VISITORS();

};
TYPEDEF_SPTR(NewtonEulerR)
#endif // NEWTONEULERRELATION_H
