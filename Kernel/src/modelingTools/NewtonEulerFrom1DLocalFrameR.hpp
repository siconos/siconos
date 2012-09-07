/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
  virtual void initComponents(Interaction& inter);

  /*Matrix converting  the absolute coordinate to the contact coordinate.*/
  SP::SimpleMatrix _Mabs_C;
  /*cross product matrices*/
  SP::SimpleMatrix _NPG1;
  SP::SimpleMatrix _NPG2;
  /*buffer matrices*/
  SP::SimpleMatrix _AUX1;
  SP::SimpleMatrix _AUX2;
private:
  void NIcomputeJachqTFromContacts(SP::NewtonEulerDS d1);
  void NIcomputeJachqTFromContacts(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2);
public:

  /** V.A. boolean _isOnCOntact ?? Why is it public members ?
  *  seems parametrize the projection algorithm
  *  the projection is done on the surface $y=0$ or on $y \geq 0$
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

  // inline void setDs1(SP::NewtonEulerDS ds1){_ds1=ds1;}
  // inline void setDs2(SP::NewtonEulerDS ds2){_ds2=ds2;}
  // inline SP::NewtonEulerDS ds1(){return _ds1;}
  // inline SP::NewtonEulerDS ds2(){return _ds2;}

  /** destructor
  */
  virtual ~NewtonEulerFrom1DLocalFrameR() {};

  virtual void computeJachq(const double time, Interaction& inter);
  virtual void computeJachqT(Interaction& inter);

  inline SP::SiconosVector pc1()
  {
    return _Pc1;
  }
  inline SP::SiconosVector pc2()
  {
    return _Pc2;
  }
  inline SP::SiconosVector nc()
  {
    return _Nc;
  }

  /** set the coordinates of first contact point
  * \param SP::SiconosVector new coordinates
  */
  void setpc1(SP::SiconosVector npc)
  {
    _Pc1 = npc;
  };

  /** set the coordinates of second contact point
  * \param SP::SiconosVector new coordinates
  */
  void setpc2(SP::SiconosVector npc)
  {
    _Pc2 = npc;
  };

  /** set the coordinates of inside normal vector at the contact point
  * \param SP::SiconosVector new coordinates
  */
  void setnc(SP::SiconosVector nnc)
  {
    _Nc = nnc;
  };

  // visitors hook
  ACCEPT_STD_VISITORS();

};
TYPEDEF_SPTR(NewtonEulerFrom1DLocalFrameR)
#endif // NEWTONEULERRIMPACT_H
