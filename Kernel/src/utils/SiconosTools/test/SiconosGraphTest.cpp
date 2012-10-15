/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "SiconosGraphTest.hpp"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosGraphTest);


void SiconosGraphTest::setUp()
{
}

void SiconosGraphTest::tearDown()
{
}

// Default constructor
void SiconosGraphTest::t1()
{

  typedef SiconosGraph < std::string, int,
          boost::no_property, boost::no_property, boost::no_property > G;

  G g;

  G::VDescriptor vd1, vd2;

  vd1 = g.add_vertex("hello");
  vd2 = g.add_vertex("goodbye");

  CPPUNIT_ASSERT(g.size() == 2);

  CPPUNIT_ASSERT(g.bundle(vd1) == "hello");
  CPPUNIT_ASSERT(g.bundle(vd2) == "goodbye");
}

void SiconosGraphTest::t2()
{

  typedef SiconosGraph < std::string, int,
          boost::no_property, boost::no_property, boost::no_property > G;

  G g;

  G::VDescriptor vd1, vd2;

  vd1 = g.add_vertex("hello");
  vd2 = g.add_vertex("goodbye");

  CPPUNIT_ASSERT(g.bundle(vd1) == "hello");
  CPPUNIT_ASSERT(g.bundle(vd2) == "goodbye");

  CPPUNIT_ASSERT(g.bundle(g.descriptor(g.bundle(vd1))) == "hello");
  CPPUNIT_ASSERT(g.bundle(g.descriptor(g.bundle(vd2))) == "goodbye");

}

void SiconosGraphTest::t3()
{

  typedef SiconosGraph < std::string, int,
          boost::no_property, boost::no_property, boost::no_property > G;

  G g;

  G::VDescriptor vd1, vd2;

  vd1 = g.add_vertex("hello");
  vd2 = g.add_vertex("goodbye");

  g.remove_vertex("hello");

  CPPUNIT_ASSERT(g.bundle(g.descriptor(g.bundle(vd2))) == "goodbye");
  g.remove_vertex("goodbye");

  CPPUNIT_ASSERT(g.size() == 0);

  vd1 = g.add_vertex("one");
  vd2 = g.add_vertex("two");

  CPPUNIT_ASSERT(g.bundle(g.descriptor(g.bundle(vd1))) == "one");
  CPPUNIT_ASSERT(g.bundle(g.descriptor(g.bundle(vd2))) == "two");

  g.remove_vertex("two");
  CPPUNIT_ASSERT(g.bundle(g.descriptor(g.bundle(vd1))) == "one");

  g.remove_vertex("one");

  CPPUNIT_ASSERT(g.size() == 0);
}

void SiconosGraphTest::t4()
{

  typedef SiconosGraph < std::string, int,
          boost::no_property, boost::no_property, boost::no_property > G;
  typedef SiconosGraph < int, std::string,
          boost::no_property, boost::no_property, boost::no_property > AG;

  G g;
  AG ag;

  G::VDescriptor vd1, vd2, vd3;

  vd1 = g.add_vertex("hello");
  vd2 = g.add_vertex("goodbye");
  vd3 = g.add_vertex("bye");

  g.add_edge(vd1, vd2, 1, ag);
  g.add_edge(vd1, vd2, 2, ag);
  g.add_edge(vd1, vd3, 3, ag);

  CPPUNIT_ASSERT(ag.size() == 3);
  CPPUNIT_ASSERT(ag.bundle(ag.descriptor(1)) == 1);
  CPPUNIT_ASSERT(ag.bundle(ag.descriptor(2)) == 2);
  CPPUNIT_ASSERT(ag.bundle(ag.descriptor(3)) == 3);
}

template<class SicGraph, class AdjointSicGraph>
struct num_inf
{
  num_inf(int n, SicGraph& sg, AdjointSicGraph& asg)
    : _n(n), _sg(sg), _asg(asg) {}
  bool operator()(typename SicGraph::EDescriptor e)
  {
    std::cout << _sg.bundle(e) << "<?" << _n << std::endl;
    if ((_sg.bundle(e) < _n) && _asg.is_vertex(_sg.bundle(e)))
    {
      CPPUNIT_ASSERT(_asg.bundle(_asg.descriptor(_sg.bundle(e))) == _sg.bundle(e)) ;

      std::cout << "removing adjoint vertex :" << _sg.bundle(e) << " <" << _n << std::endl;
      _asg.remove_vertex(_sg.bundle(e));

      CPPUNIT_ASSERT(!_asg.is_vertex(_sg.bundle(e)));
      return true;
    }
    else
    {
      return (_sg.bundle(e) < _n)
             ;
    }
  }
  int _n;
  SicGraph& _sg;
  AdjointSicGraph& _asg;
};

void SiconosGraphTest::t5()
{

  typedef SiconosGraph < std::string, int,
          boost::no_property, boost::no_property, boost::no_property > G;
  typedef SiconosGraph < int, std::string,
          boost::no_property, boost::no_property, boost::no_property > AG;

  G g;
  AG ag;

  G::VDescriptor vd1, vd2, vd3;

  vd1 = g.add_vertex("hello");
  vd2 = g.add_vertex("goodbye");
  vd3 = g.add_vertex("bye");

  std::cout << "t5:g\n";

  g.display();

  std::cout << "----\n";

  g.add_edge(vd1, vd2, 1, ag);
  g.add_edge(vd1, vd2, 2, ag);
  g.add_edge(vd1, vd3, 3, ag);

  std::cout << "t5:g+edges\n";

  g.display();

  std::cout << "----\n";

  ag.display();

  CPPUNIT_ASSERT(g.edges_number() == 3);
  CPPUNIT_ASSERT(ag.size() == 3);

  g.remove_out_edge_if(vd1, num_inf<G, AG>(3, g, ag));

  std::cout << "t5:g+remove_edge_if\n";
  g.display();
  std::cout << "----\n";
  std::cout << "t5:ag+remove_edge_if\n";
  ag.display();
  std::cout << "----\n";

  std::cout << g.edges_number() << std::endl;
  CPPUNIT_ASSERT(ag.size() == g.edges_number());

  CPPUNIT_ASSERT(g.edges_number() == 1);

  CPPUNIT_ASSERT(ag.size() == 1);

}

void SiconosGraphTest::t6()
{
  typedef SiconosGraph < std::string, int,
          boost::no_property, boost::no_property, boost::no_property > G;
  typedef SiconosGraph < int, std::string,
          boost::no_property, boost::no_property, boost::no_property > AG;

  G g;
  AG ag;

  G::VDescriptor vd1;


  vd1 = g.add_vertex("hello");

  std::cout << "t6:g\n";

  g.display();

  std::cout << "---\n";

  g.add_edge(vd1, vd1, 2, ag);
  g.add_edge(vd1, vd1, 3, ag);
  g.add_edge(vd1, vd1, 4, ag);
  g.add_edge(vd1, vd1, 5, ag);
  g.add_edge(vd1, vd1, 6, ag);
  g.add_edge(vd1, vd1, 7, ag);
  g.add_edge(vd1, vd1, 8, ag);
  g.add_edge(vd1, vd1, 9, ag);
  g.add_edge(vd1, vd1, 10, ag);
  g.add_edge(vd1, vd1, 11, ag);


  std::cout << "t6:add_edges\n";

  g.display();

  std::cout << "---\n";

  ag.display();

  std::cout << "---\n";

  CPPUNIT_ASSERT(ag.size() == g.edges_number());
  CPPUNIT_ASSERT(ag.size() == 10);

  g.remove_out_edge_if(vd1, num_inf<G, AG>(10, g, ag));

  std::cout << "t6:remove_out_edge_if < 10\n";

  g.display();

  std::cout << "---\n";

  ag.display();

  std::cout << "---\n";

  CPPUNIT_ASSERT(ag.size() == g.edges_number());

  g.remove_out_edge_if(vd1, num_inf<G, AG>(3, g, ag));

  CPPUNIT_ASSERT(ag.size() == g.edges_number());


}


void SiconosGraphTest::t7()
{

  typedef SiconosGraph < std::string, int,
          boost::no_property, boost::no_property, boost::no_property > G;
  typedef SiconosGraph < int, std::string,
          boost::no_property, boost::no_property, boost::no_property > AG;

  G g;
  AG ag;

  G::VDescriptor vd1, vd2, vd3, vd4, vd5, vd6;

  vd1 = g.add_vertex("hello");
  vd2 = g.add_vertex("goodbye");
  vd3 = g.add_vertex("bye");
  vd4 = g.add_vertex("one");
  vd5 = g.add_vertex("two");
  vd6 = g.add_vertex("three");



  g.add_edge(vd1, vd2, 1, ag);
  g.add_edge(vd2, vd3, 2, ag);
  g.add_edge(vd3, vd4, 3, ag);
  g.add_edge(vd4, vd5, 4, ag);
  g.add_edge(vd5, vd6, 5, ag);
  g.add_edge(vd1, vd5, 100, ag);
  g.add_edge(vd2, vd6, 200, ag);
  g.add_edge(vd3, vd5, 300, ag);

#ifndef NDEBUG
  CPPUNIT_ASSERT(g.state_assert());
#endif

  std::cout << "g:\n";
  g.display();

  std::cout << "ag:\n";
  ag.display();

  AG::AVIterator ui, uiend;
  std::cout << "adjacent to 100:\n";
  int tot = 0, k = 1;
  for (std11::tie(ui, uiend) = ag.adjacent_vertices(ag.descriptor(100)); ui != uiend; ++ui, k *= 10)
  {
    tot += k * ag.bundle(*ui);
  }

  CPPUNIT_ASSERT(tot == 300541);

}

void SiconosGraphTest::t8()
{
  typedef SiconosGraph < std::string, int,
          boost::no_property, boost::no_property, boost::no_property > G;
  G g;

  G::VDescriptor vd1, vd2, vd3, vd4, vd5, vd6;

  vd1 = g.add_vertex("hello");
  vd2 = g.add_vertex("goodbye");
  vd3 = g.add_vertex("bye");
  vd4 = g.add_vertex("one");
  vd5 = g.add_vertex("two");
  vd6 = g.add_vertex("three");

  CPPUNIT_ASSERT(g.bundle(vd1) == "hello");
  CPPUNIT_ASSERT(g.bundle(vd2) == "goodbye");
  CPPUNIT_ASSERT(g.bundle(vd3) == "bye");
  CPPUNIT_ASSERT(g.bundle(vd4) == "one");
  CPPUNIT_ASSERT(g.bundle(vd5) == "two");
  CPPUNIT_ASSERT(g.bundle(vd6) == "three");

}
