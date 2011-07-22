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

/*! \file SiconosGraph.hpp
  Template class to define a graph of Siconos object.

*/
#ifndef BOOST_NO_HASH
#define BOOST_NO_HASH
#endif

#ifndef SiconosGraph_H
#define SiconosGraph_H

#include <boost/config.hpp>
#include <boost/version.hpp>

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>

#if (BOOST_VERSION >= 104000)
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <boost/static_assert.hpp>

#include "SiconosSerialization.hpp"

using namespace boost;

enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };
enum graph_properties_t { graph_properties };

namespace boost
{
BOOST_INSTALL_PROPERTY(vertex, properties);
BOOST_INSTALL_PROPERTY(edge, properties);
BOOST_INSTALL_PROPERTY(graph, properties);
}



template <class V, class E, class VProperties, class EProperties, class GProperties>
class SiconosGraph
{
public:

  /* note : OutEdgeList as multisetS => cannot compile remove_out_edge_if :
     /usr/include/boost/graph/detail/adjacency_list.hpp:440: error: passing 'const ... */

  typedef adjacency_list <
  listS, listS, undirectedS,
         property < vertex_bundle_t, V,
         property < vertex_color_t ,
         default_color_type ,
         property < vertex_index_t, size_t,
         property< vertex_properties_t , VProperties > > > > ,
         property < edge_bundle_t, E,
         property < edge_color_t ,
         default_color_type ,
         property < edge_index_t, size_t,
         property< edge_properties_t , EProperties > > > > ,
         property < graph_properties_t, GProperties > >
         graph_t;

  typedef V vertex_t;

  typedef E edge_t;

  typedef typename
  graph_traits<graph_t>::edge_iterator EIterator;

  typedef typename
  graph_traits<graph_t>::vertex_iterator VIterator;

  typedef typename
  graph_traits<graph_t>::edge_descriptor EDescriptor;

  typedef typename
  graph_traits<graph_t>::vertex_descriptor VDescriptor;

  typedef typename
  graph_traits<graph_t>::out_edge_iterator OEIterator;

  typedef typename
  graph_traits<graph_t>::adjacency_iterator AVIterator;

  typedef typename
  property_map<graph_t, edge_bundle_t >::type EBundleAccess;

  typedef typename
  property_map<graph_t, vertex_bundle_t >::type VBundleAccess;

  typedef typename
  property_map<graph_t, edge_color_t >::type EColorAccess;

  typedef typename
  property_map<graph_t, vertex_color_t >::type VColorAccess;

  typedef typename
  property_map<graph_t, edge_index_t >::type EIndexAccess;

  typedef typename
  property_map<graph_t, vertex_index_t >::type VIndexAccess;

  typedef typename
  property_map<graph_t, edge_properties_t >::type EPropertiesAccess;

  typedef typename
  property_map<graph_t, vertex_properties_t >::type VPropertiesAccess;

  typedef typename
  property_map<graph_t, graph_properties_t >::type GraphPropertiesAccess;


  typedef typename std::map<V, VDescriptor> VMap;

  VMap vertex_descriptor;

protected:
  /** serialization hooks
  */
  typedef void serializable;
  template<typename Archive>
  friend void save(Archive&, SiconosGraph<V, E, VProperties, EProperties, GProperties>&, const unsigned int);
  template<typename Archive>
  friend void load(Archive&, SiconosGraph<V, E, VProperties, EProperties, GProperties>&, const unsigned int);

  graph_t g;

private:

  SiconosGraph(const SiconosGraph&);

public:

  /** default constructor
   */
  SiconosGraph()
  {
  };

  ~SiconosGraph()
  {
    g.clear();
  };


  std::pair<EDescriptor, bool>
  edge(VDescriptor u, VDescriptor v)
  {
    return boost::edge(u, v, g);
  }

  bool edge_exists(const VDescriptor& vd1, const VDescriptor& vd2)
  {
    bool ret = false;
    EDescriptor tmped;
    boost::tie(tmped, ret) = edge(vd1, vd2);

#ifndef NDEBUG
    bool check_ret = false;
    AVIterator avi, aviend;
    for (boost::tie(avi, aviend) = adjacent_vertices(vd1);
         avi != aviend; ++avi)
    {
      if (*avi == vd2)
      {
        check_ret = true;
        break;
      }
      assert(is_vertex(bundle(*avi)));
      assert(bundle(descriptor(bundle(*avi))) == bundle(*avi));
    }
    assert(ret == check_ret);
#endif

    return ret;
  }

  /* parrallel edges : edge_range needs multisetS as
     OutEdgesList and with multisetS remove_out_edges_if cannot
     compile as with listS.
     This is only needed for AdjointGraph where only 2 edges may be in
     common which correspond to the source and target in primal graph
   */
  std::pair<EDescriptor, EDescriptor>
  edges(VDescriptor u, VDescriptor v)
  {
    //    BOOST_STATIC_ASSERT((GProperties::is_adjoint_graph));

    OEIterator oei, oeiend;
    bool ifirst = false;
    bool isecond = false;
    EDescriptor first, second;
    for (boost::tie(oei, oeiend) = out_edges(u); oei != oeiend; ++oei)
    {
      if (target(*oei) == v)
      {
        if (!ifirst)
        {
          ifirst = true;
          first = *oei;
        }
        else
        {
          isecond = true;
          second = *oei;
          break;
        }
      }
    }

    if (ifirst && isecond)
    {
      if (index(first) < index(second))
      {
        return std::pair<EDescriptor, EDescriptor>(first, second);
      }
      else
      {
        return std::pair<EDescriptor, EDescriptor>(second, first);
      }
    }
    else if (ifirst)
    {
      return std::pair<EDescriptor, EDescriptor>(first, first);
    }
    else
    {
      throw(1);
    }

  }

  bool is_edge(const VDescriptor& vd1, const VDescriptor& vd2,
               const E& e_bundle)
  {
    bool found = false;
    OEIterator oei, oeiend;
    for (boost::tie(oei, oeiend) = out_edges(vd1);
         oei != oeiend; ++oei)
    {
      if (target(*oei) == vd2 && bundle(*oei) == e_bundle)
      {
        found = true;
        break;
      }
    }
    return found;
  }

  bool adjacent_vertex_exists(const VDescriptor& vd)
  {
    bool ret = false;
    VIterator vi, viend;
    for (boost::tie(vi, viend) = vertices(); vi != viend; ++vi)
    {
      assert(is_vertex(bundle(*vi)));
      assert(bundle(descriptor(bundle(*vi))) == bundle(*vi));

      ret = edge_exists(vd, *vi);
      if (ret) break;
    }
    return ret;
  }


  size_t size() const
  {
    return num_vertices(g);
  };

  size_t vertices_number() const
  {
    return num_vertices(g);
  };

  size_t edges_number() const
  {
    return num_edges(g);
  };

  inline V& bundle(const VDescriptor& vd)
  {
    return get(vertex_bundle, g)[vd];
  };

  inline E& bundle(const EDescriptor& ed)
  {
    return get(edge_bundle, g)[ed];
  };

  inline default_color_type& color(const VDescriptor& vd)
  {
    return get(vertex_color, g)[vd];
  };

  inline default_color_type& color(const EDescriptor& ed)
  {
    return get(edge_color, g)[ed];
  };

  inline GProperties& properties()
  {
    return get_property(g, graph_properties);
  };

  inline size_t& index(const VDescriptor& vd)
  {
    return get(vertex_index, g)[vd];
  };

  inline size_t& index(const EDescriptor& ed)
  {
    return get(edge_index, g)[ed];
  };

  inline VProperties& properties(const VDescriptor& vd)
  {
    return get(vertex_properties, g)[vd];
  };

  inline EProperties& properties(const EDescriptor& ed)
  {
    return get(edge_properties, g)[ed];
  };

  inline bool is_vertex(const V& vertex)
  {
    return (vertex_descriptor.find(vertex) != vertex_descriptor.end());
  }

  inline const VDescriptor& descriptor(const V& vertex)
  {
    assert(size() == vertex_descriptor.size());
    assert(vertex_descriptor.find(vertex) != vertex_descriptor.end());
    return vertex_descriptor[vertex];
  }

  inline std::pair<VIterator, VIterator> vertices()
  {
    return boost::vertices(g);
  };

  inline VIterator begin()
  {
    VIterator vi, viend;
    boost::tie(vi, viend) = vertices();
    return vi;
  }

  inline VIterator end()
  {
    VIterator vi, viend;
    boost::tie(vi, viend) = vertices();
    return viend;
  }

  inline std::pair<AVIterator, AVIterator> adjacent_vertices(const VDescriptor& vd)
  {
    return boost::adjacent_vertices(vd, g);
  };

  inline std::pair<EIterator, EIterator> edges()
  {
    return boost::edges(g);
  };

  inline std::pair<OEIterator, OEIterator> out_edges(const VDescriptor& vd)
  {
    return boost::out_edges(vd, g);
  };

  inline VDescriptor target(const EDescriptor& ed)
  {
    return boost::target(ed, g);
  };

  inline VDescriptor source(const EDescriptor& ed)
  {
    return boost::source(ed, g);
  };


  VDescriptor add_vertex(const V& vertex_bundle)
  {
    assert(vertex_descriptor.size() == size()) ;
    //    assert ( vertex_descriptor.find(vertex_bundle) == vertex_descriptor.end() );

    VDescriptor new_vertex_descriptor;
    typename VMap::iterator current_vertex_iterator =
      vertex_descriptor.find(vertex_bundle);

    if (current_vertex_iterator == vertex_descriptor.end())
    {
      new_vertex_descriptor = boost::add_vertex(g);

      assert(vertex(size() - 1, g) == new_vertex_descriptor);
      assert(size() == vertex_descriptor.size() + 1);

      vertex_descriptor[vertex_bundle] = new_vertex_descriptor;
      assert(size() == vertex_descriptor.size());

      bundle(new_vertex_descriptor) = vertex_bundle;

      assert(descriptor(vertex_bundle) == new_vertex_descriptor);
      assert(bundle(descriptor(vertex_bundle)) == vertex_bundle);

      assert(vertex_descriptor.size() == size());
      return new_vertex_descriptor;
    }
    else
    {
      assert(descriptor(vertex_bundle) == current_vertex_iterator->second);
      assert(bundle(descriptor(vertex_bundle)) == vertex_bundle);
      return current_vertex_iterator->second;
    }


    assert(false) ;
  }

  template<class G> void copy_vertex(const V& vertex_bundle, G& og)
  {

    // is G similar ?
    BOOST_STATIC_ASSERT((boost::is_same
                         <typename G::vertex_t, vertex_t>::value));
    BOOST_STATIC_ASSERT((boost::is_same
                         <typename G::edge_t, edge_t>::value));

    assert(og.is_vertex(vertex_bundle));

    VDescriptor descr = add_vertex(vertex_bundle);
    properties(descr) = og.properties(og.descriptor(vertex_bundle));

    assert(bundle(descr) == vertex_bundle);

    // edges copy as in boost::subgraph
    typename G::OEIterator ogoei, ogoeiend;
    for (boost::tie(ogoei, ogoeiend) =
           og.out_edges(og.descriptor(vertex_bundle));
         ogoei != ogoeiend; ++ogoei)
    {
      typename G::VDescriptor ognext_descr = og.target(*ogoei);

      assert(og.is_vertex(og.bundle(ognext_descr)));

      // target in graph ?
      if (is_vertex(og.bundle(ognext_descr)))
      {
        assert(bundle(descriptor(og.bundle(ognext_descr)))
               == og.bundle(ognext_descr));

        EDescriptor edescr =
          add_edge(descr, descriptor(og.bundle(ognext_descr)),
                   og.bundle(*ogoei));

        assert(bundle(edescr) == og.bundle(*ogoei));
      }
    }
  }

  void remove_vertex(const V& vertex_bundle)
  {
    assert(is_vertex(vertex_bundle));
    assert(vertex_descriptor.size() == size());
    assert(bundle(descriptor(vertex_bundle)) == vertex_bundle);


    VDescriptor vd = descriptor(vertex_bundle);

    assert(adjacent_vertices_ok());

    boost::clear_vertex(vd, g);

    assert(adjacent_vertices_ok());
    assert(!adjacent_vertex_exists(vd));

    boost::remove_vertex(vd, g);

    assert(vertex_descriptor.size() == (size() + 1));

    vertex_descriptor.erase(vertex_bundle);

    assert(adjacent_vertices_ok());
    assert(vertex_descriptor.size() == size());
    assert(!is_vertex(vertex_bundle));

    assert(state_assert());

  }

  EDescriptor add_edge(const VDescriptor& vd1,
                       const VDescriptor& vd2,
                       const E& e_bundle)
  {

    EDescriptor new_edge;
    bool inserted;

    assert(is_vertex(bundle(vd1)));
    assert(is_vertex(bundle(vd2)));

    assert(descriptor(bundle(vd1)) == vd1);
    assert(descriptor(bundle(vd2)) == vd2);

    assert(!is_edge(vd1, vd2, e_bundle));

    boost::tie(new_edge, inserted) = boost::add_edge(vd1, vd2, g);
    assert(inserted);

    bundle(new_edge) = e_bundle;

    assert(is_edge(vd1, vd2, e_bundle));

    return new_edge;
  }

  /* note : bgl self loop are a nightmare
     seg fault with self-loop ...
     https://svn.boost.org/trac/boost/ticket/4622

     try a patch r65198 on clear vertex but loop remains and stranges
     things occur...

     solution => put self loop info in VProperties (i.e outside graph
     structure)

  */


  template<class AdjointG>
  std::pair<EDescriptor, typename AdjointG::VDescriptor>
  add_edge(const VDescriptor& vd1,
           const VDescriptor& vd2,
           const E& e_bundle,
           AdjointG& ag)
  {

    // adjoint static assertions
    BOOST_STATIC_ASSERT((boost::is_same
                         <typename AdjointG::vertex_t, edge_t>::value));
    BOOST_STATIC_ASSERT((boost::is_same
                         <typename AdjointG::edge_t, vertex_t>::value));


    EDescriptor new_ed = add_edge(vd1, vd2, e_bundle);

    typename AdjointG::VDescriptor new_ve = ag.add_vertex(e_bundle);

    assert(bundle(new_ed) == ag.bundle(new_ve));
    assert(ag.size() == edges_number());

    // better to build a range [vd1,vd2] or [vd1]...
    bool endl = false;
    for (VDescriptor vdx = vd1; !endl; vdx = vd2)
    {
      assert(vdx == vd1 || vdx == vd2);

      if (vdx == vd2) endl = true;

      std::map<E, EDescriptor> Edone;

      OEIterator ied, iedend;
      for (boost::tie(ied, iedend) = out_edges(vdx);
           ied != iedend; ++ied)
      {
        if (Edone.find(bundle(*ied)) == Edone.end())
        {
          Edone[bundle(*ied)] = *ied;

          assert(source(*ied) == vdx);

          if (*ied != new_ed)
            // so this is another edge
          {
            assert(bundle(*ied) != e_bundle);

            assert(ag.bundle(ag.descriptor(bundle(*ied))) == bundle(*ied));

            assert(ag.is_vertex(bundle(*ied)));

            assert(new_ve != ag.descriptor(bundle(*ied)));

            assert(!ag.is_edge(new_ve, ag.descriptor(bundle(*ied)),
                               bundle(vdx)));

            typename AdjointG::EDescriptor aed =
              ag.add_edge(new_ve, ag.descriptor(bundle(*ied)),
                          bundle(vdx));

            assert(ag.bundle(aed) == bundle(vdx));
          }
        }

      }
    }
    assert(ag.size() == edges_number());
    return std::pair<EDescriptor, typename AdjointG::VDescriptor>(new_ed, new_ve);
  }

  void remove_edge(const EDescriptor& ed)
  {
    assert(adjacent_vertex_exists(target(ed)));
    assert(adjacent_vertex_exists(source(ed)));

    boost::remove_edge(ed, g);
    assert(state_assert());
  }

  template<class AdjointG>
  void remove_edge(const EDescriptor& ed, AdjointG& ag)
  {

    // adjoint static assertions
    BOOST_STATIC_ASSERT((boost::is_same
                         <typename AdjointG::vertex_t, edge_t>::value));
    BOOST_STATIC_ASSERT((boost::is_same
                         <typename AdjointG::edge_t, vertex_t>::value));


    assert(ag.size() == edges_number());

    assert(bundle(ed) == ag.bundle(ag.descriptor(bundle(ed))));

    remove_vertex(ag.bundle(ag.descriptor(bundle(ed))));
    remove_edge(ed);

    assert(ag.size() == edges_number());

    assert(state_assert());

  }

  /** Remove all the out-edges of vertex u for which the predicate p
   * returns true. This expression is only required when the graph
   * also models IncidenceGraph.
   */
  template<class Predicate>
  void remove_out_edge_if(const VDescriptor& vd,
                          const Predicate& pred)
  //                      Predicate pred)
  {

    BOOST_CONCEPT_ASSERT((IncidenceGraphConcept<graph_t>));
    BOOST_CONCEPT_ASSERT((MutableGraphConcept<graph_t>));

    boost::remove_out_edge_if(vd, pred, g);
    /* workaround on multisetS (tested on Disks : ok)
       multiset allows for member removal without invalidating iterators

        OEIterator oei,oeiend;
        for(boost::tie(oei,oeiend)=out_edges(vd); oei!=oeiend; ++oei)
        {
          if (pred(*oei))
          {
            remove_edge(*oei);
          }
        }
    */

    assert(state_assert());
  }

  void update_vertices_indices()
  {
    VIterator vi, viend;
    size_t i;
    for (boost::tie(vi, viend) = boost::vertices(g), i = 0;
         vi != viend; ++vi, ++i)
    {
      index(*vi) = i;
    }
  };

  void update_edges_indices()
  {
    EIterator ei, eiend;
    size_t i;
    for (boost::tie(ei, eiend) = boost::edges(g), i = 0;
         ei != eiend; ++ei, ++i)
    {
      index(*ei) = i;
    }
  };

  void update_vertices_indices(const std::vector<size_t>& indices)
  {
    VIterator vi, viend;
    size_t i;
    for (boost::tie(vi, viend) = boost::vertices(g), i = 0;
         vi != viend; ++vi, ++i)
    {
      index(*vi) = indices[i];
    }
  };

  void update_edges_indices(const std::vector<size_t>& indices)
  {
    EIterator ei, eiend;
    size_t i;
    for (boost::tie(ei, eiend) = boost::edges(g), i = 0;
         ei != eiend; ++ei, ++i)
    {
      index(*ei) = indices[i];
    }
  };


  void clear()
  {
    g.clear();
  };

  VMap vertex_descriptor_map()
  {
    return vertex_descriptor;
  };

  void display()
  {

    VIterator vi, viend;
    for (boost::tie(vi, viend) = vertices();
         vi != viend; ++vi)
    {
      std::cout << "bundle : "
                << bundle(*vi)
                << ", index : "
                << index(*vi)
                << ", color : "
                << color(*vi);
      OEIterator oei, oeiend, next;
      for (boost::tie(oei, oeiend) = out_edges(*vi);
           oei != oeiend; ++oei)
      {
        std::cout << "---"
                  << bundle(*oei)
                  << "-->"
                  << "bundle : "
                  << bundle(target(*oei))
                  << ", index : "
                  << index(target(*oei))
                  << ", color : "
                  << color(target(*oei));
      }
      std::cout << std::endl;
    }
  }

  /* debug */
#ifndef NDEBUG
  bool state_assert()
  {
    VIterator vi, viend;
    for (boost::tie(vi, viend) = vertices(); vi != viend; ++vi)
    {
      assert(is_vertex(bundle(*vi)));
      assert(bundle(descriptor(bundle(*vi))) == bundle(*vi));

      OEIterator ei, eiend;
      for (boost::tie(ei, eiend) = out_edges(*vi);
           ei != eiend; ++ei)
      {
        assert(is_vertex(bundle(target(*ei))));
        assert(source(*ei) == *vi);
      }
      AVIterator avi, aviend;
      for (boost::tie(avi, aviend) = adjacent_vertices(*vi);
           avi != aviend; ++avi)
      {
        assert(is_vertex(bundle(*avi)));
        assert(bundle(descriptor(bundle(*avi))) == bundle(*avi));
      }
    }
    return true;

  }

  bool adjacent_vertices_ok()
  {

    VIterator vi, viend;
    for (boost::tie(vi, viend) = vertices(); vi != viend; ++vi)
    {
      assert(is_vertex(bundle(*vi)));
      assert(bundle(descriptor(bundle(*vi))) == bundle(*vi));

      AVIterator avi, aviend;
      for (boost::tie(avi, aviend) = adjacent_vertices(*vi);
           avi != aviend; ++avi)
      {
        assert(is_vertex(bundle(*avi)));
        assert(bundle(descriptor(bundle(*avi))) == bundle(*avi));
      }
    }
    return true;
  }
#endif


};


#endif


