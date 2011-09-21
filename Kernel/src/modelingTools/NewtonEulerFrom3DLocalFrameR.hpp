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
#ifndef NEWTONEULERRELATIONFC3D_H
#define NEWTONEULERRELATIONFC3D_H

#include "NewtonEulerFrom1DLocalFrameR.hpp"
/** NewtonEulerFrom3DLocalFrameR
 *
 * \author O. Bonnefon
 *  \version 3.0.0.
 *  \date Dec, 2010
 *
 * This class is an interface for relation with impact and FC3D.
 * From NewtonEulerFrom1DLocalFrameR, it inherits to the computation of the jacoboian, this operator is use for the predictor of activation and deactivation of the UR.
 * The OSNSP is build using the matrix jachqT, that is computed from the point if contact pc1, pc2 and Nc.
 * Use this class consists in overload the method computeh, and children class has to set the menber pc1, pc2 and nc.
 *
 *
 */

class NewtonEulerFrom3DLocalFrameR : public NewtonEulerFrom1DLocalFrameR
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEulerFrom3DLocalFrameR);

  void FC3DcomputeJachqTFromContacts(SP::NewtonEulerDS d1);
  void FC3DcomputeJachqTFromContacts(SP::NewtonEulerDS d1, SP::NewtonEulerDS d2);
public:
  NewtonEulerFrom3DLocalFrameR(): NewtonEulerFrom1DLocalFrameR() {}

  /** destructor
   */
  virtual ~NewtonEulerFrom3DLocalFrameR() {};
  /** initialize components specific to derived classes.
    */
  virtual void initComponents();

  /*default implementation consists in multiplying jachq and T*/
  virtual void computeJachqT();

  ACCEPT_STD_VISITORS();
};
TYPEDEF_SPTR(NewtonEulerFrom3DLocalFrameR);
#endif // NEWTONEULERRELATIONFC3D_H
