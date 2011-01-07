/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
/*! \file NewtonEulerR.h

*/
#ifndef NEWTONEULERRELATIONFC3D_H
#define NEWTONEULERRELATIONFC3D_H

#include "NewtonEulerRImpact.hpp"
/** NewtonEulerRFC3D
 *
 * \author O. Bonnefon
 *  \version 3.0.0.
 *  \date Dec, 2010
 *
 * This class is an interface for relation with impact and FC3D.
 * From NewtonEulerRImpact, it inherits to the computation of the jacoboian, this operator is use for the predictor of activation and deactivation of the UR.
 * The OSNSP is build using the matrix jachqT, that is computed from the point if contact pc1, pc2 and Nc.
 * Use this class consists in overload the method computeh, and children class has to set the menber pc1, pc2 and nc.
 *
 *
 */

class NewtonEulerRFC3D : public NewtonEulerRImpact
{

protected:
public:
  NewtonEulerRFC3D(): NewtonEulerRImpact() {}

  /** destructor
   */
  virtual ~NewtonEulerRFC3D() {};
  /** initialize components specific to derived classes.
    */
  virtual void initComponents();

  /*default implementation consists in multiplying jachq and T*/
  virtual void computeJachqT();

  ACCEPT_STD_VISITORS();
};
TYPEDEF_SPTR(NewtonEulerRFC3D);
#endif // NEWTONEULERRELATIONFC3D_H
