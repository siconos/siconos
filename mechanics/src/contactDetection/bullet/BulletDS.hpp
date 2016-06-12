/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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


  /** compute and set a relative transformation of a collision object
   * \param cobj the collision object
   * \param base_translation the reference frame translation
   * \param base_rotation the reference frame rotation
   * \param offset_translation the translation from the reference frame
   * \param offset_rotation the rotation from the reference frame
   */
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
