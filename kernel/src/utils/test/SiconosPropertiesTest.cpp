// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SiconosPropertiesTest.hpp"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosPropertiesTest);

using G = siconos::graphs::SiconosGraph<std::string, std::string, boost::no_property,
                                        boost::no_property, boost::no_property>;

void SiconosPropertiesTest::setUp() {}

void SiconosPropertiesTest::tearDown() {}

/* TypeOf */
void SiconosPropertiesTest::t1() {
  auto g = std::make_shared<G>();

  auto mv = siconos::graphs::vertexProperties<int>(*g);

  auto v1 = g->add_vertex("A");
  auto v2 = g->add_vertex("B");
  auto v3 = g->add_vertex("C");
  auto v4 = g->add_vertex("D");
  auto v5 = g->add_vertex("E");

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

void SiconosPropertiesTest::t2() {
  auto g = std::make_shared<G>();

  auto mv = siconos::graphs::vertexProperties<int>(*g);

  auto v1 = g->add_vertex("A");
  auto v2 = g->add_vertex("B");
  auto v3 = g->add_vertex("C");
  auto v4 = g->add_vertex("D");
  auto v5 = g->add_vertex("E");

  g->update_vertices_indices();

  mv[v1] = 1;
  mv[v2] = 2;
  mv[v3] = 3;
  mv[v4] = 4;
  mv[v5] = 5;

  g->remove_vertex("C");
  auto vnc = g->add_vertex("new C");
  g->update_vertices_indices();

  mv[vnc] = 33;

  CPPUNIT_ASSERT(mv[v1] == 1);
  CPPUNIT_ASSERT(mv[v2] == 2);
  CPPUNIT_ASSERT(mv[v4] == 4);
  CPPUNIT_ASSERT(mv[v5] == 5);
  CPPUNIT_ASSERT(mv[vnc] == 33);
}

void SiconosPropertiesTest::t3() {
  auto g = std::make_shared<G>();

  auto mv = siconos::graphs::vertexProperties<int>(*g);

  auto v1 = g->add_vertex("A");
  auto v2 = g->add_vertex("B");
  auto v3 = g->add_vertex("C");
  auto v4 = g->add_vertex("D");
  auto v5 = g->add_vertex("E");

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

  auto vf = g->add_vertex("F");
  g->update_vertices_indices();

  mv[vf] = 333;

  CPPUNIT_ASSERT(mv[v5] == 5);
  CPPUNIT_ASSERT(mv[vf] == 333);
}
