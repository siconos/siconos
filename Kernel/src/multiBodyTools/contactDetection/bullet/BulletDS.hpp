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

#include "NewtonEulerDS.hpp"
#include "SiconosPointers.hpp"

#include "BulletSiconos.hpp"


#ifndef BulletDS_hpp
#define BulletDS_hpp

class BulletDS : public NewtonEulerDS, public boost::enable_shared_from_this<BulletDS>
{

private:
  SP::btCollisionShape _collisionShape;
  SP::btCollisionObject _collisionObject;

public:

  /** Constructor from bullet shapes enum type and generic shapeParam:
   *  collisionShape is self allocated
   */
  BulletDS(const BroadphaseNativeTypes& shape_type,
           const SP::SimpleVector& shapeParams,
           SP::SiconosVector position,
           SP::SiconosVector velocity,
           const double& mass);

  /** Constructor from collisionShape. This one should be prefered for
   *  many objects sharing same shapes
   */
  BulletDS(const SP::btCollisionShape& shape,
           SP::SiconosVector position,
           SP::SiconosVector velocity,
           const double& mass);

  /** get the collision object
  **/
  SP::btCollisionObject collisionObject() const
  {
    return _collisionObject;
  };

  /** get a shared_ptr from this
   */
  boost::shared_ptr<BulletDS> shared_ptr()
  {
    return shared_from_this();
  };

  /** visitor hook
   */
  ACCEPT_STD_VISITORS();
};

struct ForCollisionObject : public Question<SP::btCollisionObject>
{
  ANSWER(BulletDS, collisionObject());
};

TYPEDEF_SPTR(BulletDS);

#endif

