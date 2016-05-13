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
#include "SiconosPropertiesTest.hpp"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosPropertiesTest);


void SiconosPropertiesTest::setUp()
{
}

void SiconosPropertiesTest::tearDown()
{
}

/* TypeOf */
void SiconosPropertiesTest::t1()
{

  typedef SiconosGraph < std::string, std::string,
          boost::no_property, boost::no_property, boost::no_property > G;

  std11::shared_ptr<G> g(new G());


  Siconos::VertexProperties<int, G> mv = Siconos::vertexProperties<int>(*g);

  G::VDescriptor v1 = g->add_vertex("A");
  G::VDescriptor v2 = g->add_vertex("B");
  G::VDescriptor v3 = g->add_vertex("C");
  G::VDescriptor v4 = g->add_vertex("D");
  G::VDescriptor v5 = g->add_vertex("E");

  g->update_vertices_indices();

  mv[v1] = 1;
  mv[v2] = 2;
  mv[v3] = 3;
  mv[v4] = 4;
  mv[v5] = 5;

  g->remove_vertex("C");
  g->update_vertices_indices();

  CPPUNIT_ASSERT(mv[v1] == 1);
  CPPUNIT_ASSERT(mv[v2] == 2);
  CPPUNIT_ASSERT(mv[v4] == 4);
  CPPUNIT_ASSERT(mv[v5] == 5);


}

void SiconosPropertiesTest::t2()
{

  typedef SiconosGraph < std::string, std::string,
          boost::no_property, boost::no_property, boost::no_property > G;

  std11::shared_ptr<G> g(new G());


  Siconos::VertexProperties<int, G> mv = Siconos::vertexProperties<int>(*g);

  G::VDescriptor v1 = g->add_vertex("A");
  G::VDescriptor v2 = g->add_vertex("B");
  G::VDescriptor v3 = g->add_vertex("C");
  G::VDescriptor v4 = g->add_vertex("D");
  G::VDescriptor v5 = g->add_vertex("E");

  g->update_vertices_indices();

  mv[v1] = 1;
  mv[v2] = 2;
  mv[v3] = 3;
  mv[v4] = 4;
  mv[v5] = 5;

  g->remove_vertex("C");
  G::VDescriptor vnc = g->add_vertex("new C");
  g->update_vertices_indices();

  mv[vnc] = 33;

  CPPUNIT_ASSERT(mv[v1] == 1);
  CPPUNIT_ASSERT(mv[v2] == 2);
  CPPUNIT_ASSERT(mv[v4] == 4);
  CPPUNIT_ASSERT(mv[v5] == 5);
  CPPUNIT_ASSERT(mv[vnc] == 33);

}

void SiconosPropertiesTest::t3()
{
  typedef SiconosGraph < std::string, std::string,
          boost::no_property, boost::no_property, boost::no_property > G;

  std11::shared_ptr<G> g(new G());


  Siconos::VertexProperties<int, G> mv = Siconos::vertexProperties<int>(*g);

  G::VDescriptor v1 = g->add_vertex("A");
  G::VDescriptor v2 = g->add_vertex("B");
  G::VDescriptor v3 = g->add_vertex("C");
  G::VDescriptor v4 = g->add_vertex("D");
  G::VDescriptor v5 = g->add_vertex("E");

  g->update_vertices_indices();

  mv[v1] = 1;
  mv[v2] = 2;
  mv[v3] = 3;
  mv[v4] = 4;
  mv[v5] = 5;

  g->remove_vertex("A");
  g->remove_vertex("B");
  g->remove_vertex("C");
  g->remove_vertex("D");

  G::VDescriptor vf = g->add_vertex("F");
  g->update_vertices_indices();

  mv[vf] = 333;

  CPPUNIT_ASSERT(mv[v5] == 5);
  CPPUNIT_ASSERT(mv[vf] == 333);

}
