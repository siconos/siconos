// -*- c++ -*-
// SWIG interface for Siconos Mechanics/ContactDetection/Bullet
%module(package="siconos.mechanics.collision", directors="1", allprotected="1") bullet

 // serialization not yet implemented for bullet
#undef WITH_IO
#undef WITH_SERIALIZATION

#ifdef BT_USE_DOUBLE_PRECISION
%{
#define BT_USE_DOUBLE_PRECISION 1
%}
#endif

%include native.i
%include base.i

// due to undefined private copy constructors
%feature("notabstract") BulletTimeStepping;

// do not wrap visitor visit : this lead to a huge amount of wrapper
// code generation and this fail at compile time on shared_ptr freearg
%ignore visit;

// ignores mostly because not defined in <name>.h
%shared_ptr(btCapsuleShapeX);
%ignore btCapsuleShapeZ;
%ignore btConeShapeX;
%ignore btConeShapeZ;
%shared_ptr(btCylinderShapeX);
%shared_ptr(btCylinderShapeZ);
%ignore btBU_Simplex1to4;
%ignore m_vertices1;

%ignore btVector4;

%ignore btVector3::m_floats;
%ignore btFace::m_plane;

#undef PY_REGISTER_BULLET_COLLISION_DETECTION
%define PY_REGISTER_BULLET_COLLISION_DETECTION(X)
%inline
%{
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif
#include <BulletCollision/CollisionShapes/X.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic pop
#endif
%}
%shared_ptr(X);
%include "BulletCollision/CollisionShapes/X.h";
%enddef

#undef PY_REGISTER_BULLET_NARROW_PHASE_COLLISION_DETECTION
%define PY_REGISTER_BULLET_NARROW_PHASE_COLLISION_DETECTION(X)
%inline
%{
#include <BulletCollision/NarrowPhaseCollision/X.h>
%}
%shared_ptr(X);
%include "BulletCollision/NarrowPhaseCollision/X.h";
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
%import LinearMath/btScalar.h

PY_REGISTER_BULLET_LINEAR_MATH(btVector3);
PY_REGISTER_BULLET_LINEAR_MATH(btQuadWord);
PY_REGISTER_BULLET_LINEAR_MATH(btQuaternion);
PY_REGISTER_BULLET_LINEAR_MATH(btMatrix3x3);
PY_REGISTER_BULLET_LINEAR_MATH(btTransform);

%{
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
%}
%shared_ptr(btDefaultCollisionConfiguration);
%include "BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h"

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
#include <BulletCollision/BroadphaseCollision/btBroadphaseInterface.h>
%}
%shared_ptr(btBroadphaseInterface);
%include "BulletCollision/BroadphaseCollision/btBroadphaseInterface.h"

%{
#include <BulletCollision/CollisionDispatch/btCollisionWorld.h>
%}
%shared_ptr(btCollisionWorld);
%include "BulletCollision/CollisionDispatch/btCollisionWorld.h"


%shared_ptr(std::vector< std11::shared_ptr<btCollisionObject> >);
%template (collisionObjects) std::vector< std11::shared_ptr< btCollisionObject > >;

//%shared_ptr(std::vector< std11::shared_ptr<btCollisionShape> >);

PY_REGISTER_BULLET_NARROW_PHASE_COLLISION_DETECTION(btManifoldPoint);
PY_REGISTER_BULLET_NARROW_PHASE_COLLISION_DETECTION(btPersistentManifold);

// For BulletR::BulletR()
REF_PTR(btManifoldPoint)

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
PY_REGISTER_BULLET_COLLISION_DETECTION(btStridingMeshInterface);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleIndexVertexArray);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleIndexVertexMaterialArray);
PY_REGISTER_BULLET_COLLISION_DETECTION(btTriangleInfoMap);
PY_REGISTER_BULLET_COLLISION_DETECTION(btUniformScalingShape);



%{
#include <BulletCollision/Gimpact/btGImpactShape.h>
%}
%shared_ptr(btTetrahedronShapeEx);
%shared_ptr(btGImpactShapeInterface);
%shared_ptr(btGImpactCompoundShape);
%shared_ptr(btGImpactMeshShapePart);
%shared_ptr(btGImpactShape);
%shared_ptr(btGImpactMeshShape);
%include "BulletCollision/Gimpact/btGImpactShape.h"

// force the definition of SWIGTYPE_p_Interaction...
typedef Interaction Interaction;

%include "BulletSiconosFwd.hpp"
PY_FULL_REGISTER(BulletR);
PY_FULL_REGISTER(BulletDS);
PY_FULL_REGISTER(BulletSpaceFilter);
PY_FULL_REGISTER(BulletTimeStepping);
PY_FULL_REGISTER(BulletTimeSteppingDirectProjection);
PY_FULL_REGISTER(BulletWeightedShape);
PY_FULL_REGISTER(BulletFrom1DLocalFrameR);




%inline
{
  SP::BulletDS cast_BulletDS(SP::DynamicalSystem ds)
  {
    return std11::dynamic_pointer_cast<BulletDS>(ds);
  };

  SP::BulletR cast_BulletR(SP::Relation rel)
  {
    return std11::dynamic_pointer_cast<BulletR>(rel);
  };

  extern bool gContactCalcArea3Points;


  extern  btScalar gContactBreakingThreshold;

  void set_gContactBreakingThreshold(double x)
  {
    gContactBreakingThreshold = btScalar(x);
  };
  btScalar get_gContactBreakingThreshold()
  {
    return gContactBreakingThreshold;
  };

}




%extend btCollisionObject
{

  size_t __hash__()
  {
    return (size_t) $self;
  };

}

%extend btTriangleIndexVertexArray
{
  btTriangleIndexVertexArray(PyObject *o1, PyObject *o2)
  {
    int is_new_object1=0;
    int is_new_object2=0;
    PyArrayObject* points = (PyArrayObject*) o1;
    PyArrayObject* indices = (PyArrayObject*) o2;


    int num_triangles = array_size(indices,0);
    int num_vertices = array_size(points,0);

    btTriangleIndexVertexArray* index =
      new btTriangleIndexVertexArray(num_triangles, (int *) array_data(indices),
                                     3 * sizeof(int),
                                     num_vertices, (btScalar *) array_data(points),
                                     3 * sizeof(btScalar));

     // python mem management
    if(is_new_object1 && points)
    {
      Py_DECREF(points);
    }

    if(is_new_object2 && indices)
    {
      Py_DECREF(indices);
    }

    return index;
  }
}

// New Bullet stuff

%include base.i

PY_REGISTER_WITHOUT_HEADER(SiconosBulletOptions);
PY_REGISTER_WITHOUT_HEADER(SiconosBulletStatistics);
PY_FULL_REGISTER(SiconosBulletCollisionManager);
