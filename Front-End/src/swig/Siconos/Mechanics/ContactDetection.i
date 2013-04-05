%module(directors="1", allprotected="1") ContactDetection

%include start.i

%include path.i

%include handleException.i

%include sharedPointers.i

%include KernelTypes.i
%{
#include <SiconosKernel.hpp>
%}

%import Kernel.i

%include pyRegister.i

%fragment("NumPy_Fragments");
PY_REGISTER(SpaceFilter);                                             

PY_REGISTER(SiconosBodies);
// yes, undefined private copy constructors
%feature("notabstract") BulletTimeStepping;


#ifdef WITH_BULLET

// do not wrap visitor visit : this lead to a huge amount of wrapper
// code generation and this fail at compile time on shared_ptr freearg
%ignore visit;

// ignores mostly because not defined in <name>.h
%shared_ptr(btCapsuleShapeX);
%ignore btCapsuleShapeZ;
%ignore btConeShapeX;
%ignore btConeShapeZ;
%ignore btCylinderShapeX;
%ignore btCylinderShapeZ;
%ignore btBU_Simplex1to4;
%ignore m_vertices1;

%ignore btVector4;

%ignore btVector3::m_floats;
%ignore btFace::m_plane;

#undef PY_REGISTER_BULLET_COLLISION_DETECTION
%define PY_REGISTER_BULLET_COLLISION_DETECTION(X)
%inline
%{
#include <BulletCollision/CollisionShapes/X.h>
%}
%shared_ptr(X);
%include "BulletCollision/CollisionShapes/X.h";
%enddef

#undef PY_REGISTER_BULLET_LINEAR_MATH
%define PY_REGISTER_BULLET_LINEAR_MATH(X)
%inline
%{
#include <LinearMath/X.h>
%}
%shared_ptr(X);
%include "LinearMath/X.h";
%enddef

%shared_ptr(btCollisionShape);
%shared_ptr(btConvexShape);
%shared_ptr(btConvexInternalShape);
%shared_ptr(btConvexInternalAabbCachingShape);
%shared_ptr(btPolyhedralConvexShape);
%shared_ptr(btPolyhedralConvexAabbCachingShape);
%shared_ptr(btConvexHullShape);


%{
#include <LinearMath/btScalar.h>
%}
%include LinearMath/btScalar.h


PY_REGISTER_BULLET_LINEAR_MATH(btVector3);
PY_REGISTER_BULLET_LINEAR_MATH(btMatrix3x3);
PY_REGISTER_BULLET_LINEAR_MATH(btTransform);

%{
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>
%}
%shared_ptr(btCollisionObject);
%include "BulletCollision/CollisionDispatch/btCollisionObject.h"

%{
#include <BulletCollision/BroadphaseCollision/btDispatcher.h>
%}
%shared_ptr(btDispatcher);
%include "BulletCollision/BroadphaseCollision/btDispatcher.h"

%{
#include <BulletCollision/CollisionDispatch/btCollisionWorld.h>
%}
%shared_ptr(btCollisionWorld);
%include "BulletCollision/CollisionDispatch/btCollisionWorld.h"



PY_REGISTER_BULLET_COLLISION_DETECTION(btCollisionShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexInternalShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvex2dShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btPolyhedralConvexShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexHullShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexPointCloudShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexPolyhedron);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexTriangleMeshShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConcaveShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btEmptyShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btCompoundShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleMesh);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleMeshShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btBox2dShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btBoxShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btCapsuleShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConeShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btCylinderShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btHeightfieldTerrainShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btMaterial);
PY_REGISTER_BULLET_COLLISION_DETECTION(btMinkowskiSumShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btSphereShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btMultiSphereShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btMultimaterialTriangleMeshShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btOptimizedBvh);
PY_REGISTER_BULLET_COLLISION_DETECTION(btScaledBvhTriangleMeshShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btShapeHull);
PY_REGISTER_BULLET_COLLISION_DETECTION(btStaticPlaneShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btStridingMeshInterface);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTetrahedronShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleBuffer);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleCallback);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleIndexVertexArray);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleIndexVertexMaterialArray);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleInfoMap);
PY_REGISTER_BULLET_COLLISION_DETECTION(btUniformScalingShape);

%include "BulletSiconos.hpp"
PY_REGISTER(BulletR);
PY_REGISTER(BulletDS);
PY_REGISTER(BulletSpaceFilter);
PY_REGISTER(BulletTimeStepping);
PY_REGISTER(BulletTimeSteppingProjectOnConstraints);
PY_REGISTER(BulletWeightedShape);
PY_REGISTER(BulletFrom1DLocalFrameR);

#endif
