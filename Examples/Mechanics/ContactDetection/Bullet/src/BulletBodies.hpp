#ifndef BulletBodies_hpp
#define BulletBodies_hpp
#include <SiconosBodies.hpp>

#include <SiconosKernel.hpp>
#include <bullet/LinearMath/btVector3.h>
#include <bullet/BulletCollision/CollisionShapes/btBoxShape.h>
#include <bullet/BulletCollision/CollisionShapes/btCapsuleShape.h>
#include <bullet/BulletCollision/CollisionShapes/btCylinderShape.h>
#include <bullet/BulletCollision/CollisionShapes/btSphereShape.h>
#include <bullet/BulletCollision/CollisionShapes/btConeShape.h>
#include <bullet/BulletCollision/CollisionShapes/btConvexHullShape.h>
#include <bullet/BulletCollision/CollisionShapes/btMultiSphereShape.h>
#include <bullet/BulletCollision/CollisionDispatch/btCollisionObject.h>
#include <bullet/BulletCollision/CollisionDispatch/btCollisionWorld.h>

#include <BulletDS.hpp>

class BulletBodies : public SiconosBodies
{

public:

  void init();

};
#endif
