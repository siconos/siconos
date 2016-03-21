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

public:

  /** constructor from a BulletWeightedShape
   * \param weightedShape the main bullet collision shape with an associated mass matrix
   * \param position the initial position (vector length: 3)
   * \param velocity the inital velocity  (vector length: 3)
   * \param relative_position relative position of the main bullet
   *        collision shape (default (0,0,0))
   * \param relative_orientation relative orientation (given as a
   *        quaternion) of the main buller collision shape (default
   *        (0, 1, 0, 0))
   * \param group the collision group (default 0)
   */
  BulletDS(SP::BulletWeightedShape weightedShape,
           SP::SiconosVector position,
           SP::SiconosVector velocity,
           SP::SiconosVector relative_position = SP::SiconosVector(),
           SP::SiconosVector relative_orientation = SP::SiconosVector(),
           int group=0);


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

  /** add a collision object
   * \param cobj the collision object
   * \param pos the position (x, y, z) relative to the center of mass
   * \param ori the orientation (quaternion) relative to the center of mass
   * \param group the contactor group id
   */
  void addCollisionObject(SP::btCollisionObject cobj,
                          SP::SiconosVector pos,
                          SP::SiconosVector ori,
                          int group);

  /** add a collision shape
   * \param shape the collision shape
   * \param pos the position (x, y, z) relative to the center of mass
   * \param ori the orientation (quaternion) relative to the center of mass
   * \param group the contactor group id
  */
  void addCollisionShape(SP::btCollisionShape shape,
                         SP::SiconosVector pos,
                         SP::SiconosVector ori,
                         int group);

  /** update Bullet collision objects positions and orientations
   */
  void updateCollisionObjects() const;


  static void setRelativeTransform(SP::btCollisionObject cobj,
                                   const btVector3& base_translation,
                                   const btQuaternion& base_rotation,
                                   const btVector3& offset_translation,
                                   const btQuaternion& offset_rotation);



  /** visitor hook
   */
  ACCEPT_BASE_STD_VISITORS(NewtonEulerDS);

  /** return the shared pointer associated
   * \return std11::shared_ptr<BulletDS>
   */
  std11::shared_ptr<BulletDS> shared_ptr()
  {
    return shared_from_this();
  };


private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BulletDS);

  SP::BulletWeightedShape _weightedShape;
  SP::CollisionObjects _collisionObjects;

};

struct ForCollisionObjects : public Question<SP::CollisionObjects>
{
  ANSWER(BulletDS, collisionObjects());
};


struct ForWeightedShape : public Question<SP::BulletWeightedShape>
{
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
