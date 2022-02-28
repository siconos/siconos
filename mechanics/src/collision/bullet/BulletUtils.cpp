/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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


#include "BulletUtils.hpp"

#define DEBUG_STDOUT
#define DEBUG_MESSAGES 1
#include "siconos_debug.h"

#include <BulletCollision/GImpact/btGImpactShape.h>

void display_info_collision_object(const btCollisionObject* collisionObject)
{
  const btCollisionShape* collisionShape = collisionObject->getCollisionShape();
  printf("collisionShape : %p\n", collisionShape);
  printf("collisionShape->getShapeType(): %i\n", collisionShape->getShapeType());
  printf("collisionShape->getName(): %s\n", collisionShape->getName());


  if (collisionShape->isPolyhedral())
  {
    printf("isPolyhedral() true \n");
  }
  else
    printf("isPolyhedral() false \n");

  if (collisionShape->getShapeType() == 25)
  {
    printf("GImpactMesh shape type\n");
    btGImpactMeshShape * gimpact_mesh_shape = (btGImpactMeshShape *)collisionShape;
    btStridingMeshInterface * striding_mesh = gimpact_mesh_shape->getMeshInterface();
    btTriangleIndexVertexArray * triangle_mesh = (btTriangleIndexVertexArray *) striding_mesh;
    IndexedMeshArray& mesh_array =  triangle_mesh->getIndexedMeshArray();
    for (int i =0; i < triangle_mesh->getNumSubParts(); i++)
    {
      printf(" mesh_array[%i] : number of triangles = %i\n", i, mesh_array[i].m_numTriangles);
      int * triangleIndexBase = (int*) mesh_array[i].m_triangleIndexBase;
      int k =0;
      for (int t =0; t < mesh_array[i].m_numTriangles; t++)
      {
        printf("              vertex indices  of triangle %i : %i\t %i\t %i\n", t, triangleIndexBase[k], triangleIndexBase[k+1],triangleIndexBase[k+2]);
        k=k+3;
      }
      printf("               : number of vertices = %i\n", mesh_array[i].m_numVertices);
      btScalar * vertexBase = (btScalar*) mesh_array[i].m_vertexBase;
      k =0;
      for (int v =0; v < mesh_array[i].m_numVertices; v++)
      {
        printf("              vertices  %i : %e\t %e\t %e\n", v, vertexBase[k], vertexBase[k+1], vertexBase[k+2]);
        k=k+3;
      }
    }
    printf("gimpact_mesh_shape->getMeshPartCount() = % i \n", gimpact_mesh_shape->getMeshPartCount());
    for (int mesh_part_index=0; mesh_part_index < gimpact_mesh_shape->getMeshPartCount(); mesh_part_index++ )
    {
      btGImpactMeshShapePart* mesh_part =gimpact_mesh_shape->getMeshPart(mesh_part_index);
    }
  }
  //getchar();
}

void display_info_manifold(const btPersistentManifold& manifold)
{

  printf("-------- manifold : %p\n",  &manifold);

  const btCollisionObject* body0 = manifold.getBody0();
  printf("-------- m_body0 : %p\n", body0);
  display_info_collision_object(body0);


  const btCollisionObject* body1 = manifold.getBody1();
  printf("-------- m_body1 : %p\n", body1);
  display_info_collision_object(body1);

  printf("Number of contact points (m_cachedPoints) =%i \n", manifold.getNumContacts());
  for (int index =0; index < manifold.getNumContacts(); index++)
  {
    const btManifoldPoint& point =  manifold.getContactPoint(index);
    printf("   contact point number %i : \n", index);
    btVector3  pA = point.m_localPointA;
    printf("   pA x , y, x : %e\t, %e\t, %e\t \n", pA.x(), pA.y(), pA.z());
    btVector3  pB = point.m_localPointB;
    printf("   pB x , y, x : %e\t, %e\t, %e\t \n", pB.x(), pB.y(), pB.z());
  }
}
