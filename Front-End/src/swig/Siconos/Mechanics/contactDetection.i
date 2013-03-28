%module(directors="1", allprotected="1") contactDetection

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

#ifdef WITH_BULLET
// ignores mostly because not defined in <name>.h

// do not wrap visitor visit : this lead to a huge amount of wrapper
// code generation and this fail at compile time on shared_ptr freearg
%ignore visit;

%ignore btCapsuleShapeX;
%ignore btCapsuleShapeZ;
%ignore btConeShapeX;
%ignore btConeShapeZ;
%ignore btCylinderShapeX;
%ignore btCylinderShapeZ;
%ignore btConvexInternalAabbCachingShape;
%ignore btPolyhedralConvexAabbCachingShape;
%ignore btBU_Simplex1to4;
%ignore m_vertices1;

%ignore btVector3::serialize;
%ignore btVector3::deSerialize;

%ignore btVector4;

#undef PY_REGISTER_BULLET_COLLISION_DETECTION
%define PY_REGISTER_BULLET_COLLISION_DETECTION(X)
%inline
%{
#include <BulletCollision/CollisionShapes/X.h>
%}
TYPEDEF_SPTR(X);
%include "BulletCollision/CollisionShapes/X.h";
%enddef

#undef PY_REGISTER_BULLET_LINEAR_MATH
%define PY_REGISTER_BULLET_LINEAR_MATH(X)
%inline
%{
#include <LinearMath/X.h>
%}
TYPEDEF_SPTR(X);
%include "LinearMath/X.h";
%enddef

PY_REGISTER(BulletR);
PY_REGISTER(BulletDS);
PY_REGISTER(BulletSpaceFilter);
PY_REGISTER(BulletTimeStepping);
PY_REGISTER(BulletTimeSteppingProjectOnConstraints);
PY_REGISTER(BulletWeightedShape);
PY_REGISTER(BulletFrom1DLocalFrameR);
PY_REGISTER_BULLET_LINEAR_MATH(btScalar);
PY_REGISTER_BULLET_LINEAR_MATH(btVector3);

PY_REGISTER_BULLET_COLLISION_DETECTION(btCollisionShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btCollisionMargin);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexInternalShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvex2dShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexPointCloudShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexHullShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexPolyhedron);
PY_REGISTER_BULLET_COLLISION_DETECTION(btConvexTriangleMeshShape);
PY_REGISTER_BULLET_COLLISION_DETECTION(btPolyhedralConvexShape);
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
#endif
