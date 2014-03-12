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

#ifndef BulletDS_hpp
#define BulletDS_hpp

#include "BulletSiconosFwd.hpp"
#include "NewtonEulerDS.hpp"

#include <Question.hpp>

class BulletDS : public NewtonEulerDS, public std11::enable_shared_from_this<BulletDS>
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BulletDS);

  SP::BulletWeightedShape _weightedShape;
  SP::CollisionObjects _collisionObjects;

public:

  /** Constructor from a BulletBodyShape
   */
  BulletDS(SP::BulletWeightedShape weightedShape,
           SP::SiconosVector position,
           SP::SiconosVector velocity);


  /** get the number of collision objects
      \return unsigned int the number of collision objects
   */
  unsigned int numberOfCollisionObjects() const;

  /** get collision objects
      \return SP::CollisionObjects a pointer on collision objects
  **/
  SP::CollisionObjects collisionObjects() const;

  /** get the shape responsible of the mass matrix
      \return SP::BulletWeightedShape a pointer on the shape
  **/
  SP::BulletWeightedShape weightedShape() const
  {
    return _weightedShape;
  };

  /** update Bullet collision objects positions and orientations
   */
  void updateCollisionObjects() const;

  /** visitor hook
   */
  ACCEPT_STD_VISITORS();

  /** return the shared pointer associated
   * \return std11::shared_ptr<BulletDS>
   */
  std11::shared_ptr<BulletDS> shared_ptr()
  {
    return shared_from_this();
  };

};

struct ForCollisionObjects : public Question<SP::CollisionObjects>
{
  using SiconosVisitor::visit;
  ANSWER(BulletDS, collisionObjects());
};


struct ForWeightedShape : public Question<SP::BulletWeightedShape>
{
  using SiconosVisitor::visit;

  ANSWER(BulletDS, weightedShape());
};

struct UpdateCollisionObjects : public SiconosVisitor
{
  using SiconosVisitor::visit;

  void visit(const BulletDS& bds)
  {
    bds.updateCollisionObjects();
  }
};
#endif

