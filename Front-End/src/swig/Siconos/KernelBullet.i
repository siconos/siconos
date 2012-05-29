%{
#include "BulletSiconos.hpp"
#include "BulletSpaceFilter.hpp"
#include "BulletR.hpp"
#include "BulletDS.hpp"
#include "BulletFrom1DLocalFrameR.hpp"
#include "BulletTimeStepping.hpp"
#include "BulletTimeSteppingProjectOnConstraints.hpp"
#include "BulletWeightedShape.hpp"

#include "LinearMath/btVector3.h"
#include "btBulletCollisionCommon.h"

#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btCollisionMargin.h"
#include "BulletCollision/CollisionShapes/btConvexShape.h"
#include "BulletCollision/CollisionShapes/btConvexInternalShape.h"
#include "BulletCollision/CollisionShapes/btConvex2dShape.h"
#include "BulletCollision/CollisionShapes/btConvexPointCloudShape.h"
#include "BulletCollision/CollisionShapes/btConvexHullShape.h"
#include "BulletCollision/CollisionShapes/btConvexPolyhedron.h"
#include "BulletCollision/CollisionShapes/btConvexTriangleMeshShape.h"
#include "BulletCollision/CollisionShapes/btPolyhedralConvexShape.h"

#include "BulletCollision/CollisionShapes/btConcaveShape.h"

#include "BulletCollision/CollisionShapes/btEmptyShape.h"

#include "BulletCollision/CollisionShapes/btCompoundShape.h"

#include "BulletCollision/CollisionShapes/btTriangleShape.h"

#include "BulletCollision/CollisionShapes/btTriangleMesh.h"
#include "BulletCollision/CollisionShapes/btTriangleMeshShape.h"


#include "BulletCollision/CollisionShapes/btBox2dShape.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletCollision/CollisionShapes/btCapsuleShape.h"

#include "BulletCollision/CollisionShapes/btConeShape.h"


#include "BulletCollision/CollisionShapes/btCylinderShape.h"


#include "BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h"
#include "BulletCollision/CollisionShapes/btMaterial.h"

#include "BulletCollision/CollisionShapes/btMinkowskiSumShape.h"

#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btMultiSphereShape.h"
#include "BulletCollision/CollisionShapes/btMultimaterialTriangleMeshShape.h"
#include "BulletCollision/CollisionShapes/btOptimizedBvh.h"
#include "BulletCollision/CollisionShapes/btScaledBvhTriangleMeshShape.h"
#include "BulletCollision/CollisionShapes/btShapeHull.h"

#include "BulletCollision/CollisionShapes/btStaticPlaneShape.h"
#include "BulletCollision/CollisionShapes/btStridingMeshInterface.h"
#include "BulletCollision/CollisionShapes/btTetrahedronShape.h"
#include "BulletCollision/CollisionShapes/btTriangleBuffer.h"
#include "BulletCollision/CollisionShapes/btTriangleCallback.h"
#include "BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h"
#include "BulletCollision/CollisionShapes/btTriangleIndexVertexMaterialArray.h"
#include "BulletCollision/CollisionShapes/btTriangleInfoMap.h"
#include "BulletCollision/CollisionShapes/btUniformScalingShape.h"



%}
#define SIMD_FORCE_INLINE inline

%ignore serialize;
%ignore deSerialize;
%import "LinearMath/btScalar.h"
%include "btBulletCollisionCommon.h"


