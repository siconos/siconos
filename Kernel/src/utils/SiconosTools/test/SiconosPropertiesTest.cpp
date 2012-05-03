/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

  boost::shared_ptr<G> g(new G());


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

  boost::shared_ptr<G> g(new G());


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

  boost::shared_ptr<G> g(new G());


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
