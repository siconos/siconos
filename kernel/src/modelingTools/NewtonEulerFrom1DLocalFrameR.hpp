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
/*! \file NewtonEulerFrom1DLocalFrameR.hpp

 */
#ifndef NEWTONEULERIMPACT_H
#define NEWTONEULERIMPACT_H

#include "NewtonEulerR.hpp"
#include "NewtonEulerDS.hpp"


/** NewtonEulerFrom1DLocalFrameR
 *
 * \author O. Bonnefon
 *  \version 3.0.0.
 *  \date Dec, 2010
 *
 * This class is an interface for a relation with impact.  It
 * implements the computation of the jacoboian of h from the points of
 * contacts and the normal.  Use this class consists in overloading
 * the method computeh, by setting the menber pc1, pc2, nc and y.  The
 * matrix jachq is used both for the building of the OSNSP (with T)
 * and for the predictor of activation of deactivation of the Interaction.
 *
 */


class NewtonEulerFrom1DLocalFrameR : public NewtonEulerR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEulerFrom1DLocalFrameR);

  /*Point of contact*/
  SP::SiconosVector _Pc1;
  SP::SiconosVector _Pc2;

  /*Inward Normal at the contact */
  SP::SiconosVector _Nc;

  // /*because set is not sorted!*/
  // SP::NewtonEulerDS _ds1;
  // SP::NewtonEulerDS _ds2;
  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

  /*Matrix converting  the absolute coordinate to the contact coordinate.*/
  SP::SimpleMatrix _Mabs_C;
  /* Matrix converting */
  SP::SimpleMatrix _MObjToAbs;
  /*cross product matrices*/
  SP::SimpleMatrix _NPG1;
  SP::SimpleMatrix _NPG2;
  /*buffer matrices*/
  SP::SimpleMatrix _AUX1;
  SP::SimpleMatrix _AUX2;
private:
  void NIcomputeJachqTFromContacts(SP::SiconosVector q1);
  void NIcomputeJachqTFromContacts(SP::SiconosVector q1, SP::SiconosVector q2);
public:

  /** V.A. boolean _isOnCOntact ?? Why is it public members ?
  *  seems parametrize the projection algorithm
  *  the projection is done on the surface \f$y=0\f$ or on \f$y \geq 0\f$
  */
  bool _isOnContact;

  /** constructorx
  */
  NewtonEulerFrom1DLocalFrameR():
    NewtonEulerR(), _Pc1(new SiconosVector(3)), _Pc2(new SiconosVector(3)),
    _Nc(new SiconosVector(3))
  {
    /*_ds1=NULL;_ds2=NULL;*/
  }

  /** destructor
  */
  virtual ~NewtonEulerFrom1DLocalFrameR() {};

  virtual void computeJachq(double time, Interaction& inter, VectorOfBlockVectors& DSlink);
  virtual void computeJachq(double time, Interaction& inter, SP::SiconosVector q1, SP::SiconosVector q2=SP::SiconosVector());
  
  /* Default implementation consists in multiplying jachq and T (see NewtonEulerR::computeJachqT)
   * but here we compute the operator from the the contact point locations
   * and the local frame at contact
   *  \param inter interaction that owns the relation
   *  \param DSlink the container of the link to DynamicalSystem attributes
   */
  virtual void computeJachqT(Interaction& inter, VectorOfBlockVectors& DSlink );
  
  virtual void computeJachqT(Interaction& inter, SP::SiconosVector q1, SP::SiconosVector q2);
  
  inline SP::SiconosVector pc1() const
  {
    return _Pc1;
  }
  inline SP::SiconosVector pc2() const
  {
    return _Pc2;
  }
  inline SP::SiconosVector nc() const
  {
    return _Nc;
  }

  /** set the coordinates of first contact point
  * \param npc new coordinates
  */
  void setpc1(SP::SiconosVector npc)
  {
    _Pc1 = npc;
  };

  /** set the coordinates of second contact point
  * \param npc new coordinates
  */
  void setpc2(SP::SiconosVector npc)
  {
    _Pc2 = npc;
  };

  /** set the coordinates of inside normal vector at the contact point
  * \param nnc new coordinates
  */
  void setnc(SP::SiconosVector nnc)
  {
    _Nc = nnc;
  };

  // visitors hook
  ACCEPT_STD_VISITORS();

};
#endif // NEWTONEULERRIMPACT_H
