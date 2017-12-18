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

/*! \file SiconosGraph.hpp
  Template class to define a graph of Siconos object.

  Note: this need documentation

*/

#ifndef SICONOS_GRAPH_HPP
#define SICONOS_GRAPH_HPP

#ifndef BOOST_NO_HASH
#define BOOST_NO_HASH
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif

#include <boost/config.hpp>
#include <boost/version.hpp>

#include <SiconosConfig.h>
#if defined(SICONOS_STD_UNORDERED_MAP) && !defined(SICONOS_USE_MAP_FOR_HASH)
#include <unordered_map>
#else
#include <map>
#endif

#include <limits>

/* gccxml 0.9 complains about ambiguous usage of size_t or std::size_t
 * in some boost headers, so we specify which one we want. It seems
 * that there is no difference anyway:
 * http://stackoverflow.com/questions/5813700/difference-between-size-t-and-stdsize-t */
using std::size_t;

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/directed_graph.hpp>

#if (BOOST_VERSION >= 104000)
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <boost/static_assert.hpp>

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#if defined(SICONOS_STD_TUPLE) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
namespace std11 = std;
#else
namespace std11 = boost;
#endif

#include "SiconosSerialization.hpp"

enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };
enum graph_properties_t { graph_properties };

// We need those with Boost >= 1.51
// vertex_bundle_t and edge_bundle_t are already defined in
// boost/pending/property.hpp
enum vertex_siconos_bundle_t { vertex_siconos_bundle };
enum edge_siconos_bundle_t { edge_siconos_bundle };

namespace boost
{
BOOST_INSTALL_PROPERTY(vertex, properties);
BOOST_INSTALL_PROPERTY(edge, properties);
BOOST_INSTALL_PROPERTY(graph, properties);
BOOST_INSTALL_PROPERTY(vertex, siconos_bundle);
BOOST_INSTALL_PROPERTY(edge, siconos_bundle);
}



template < class V, class E, class VProperties,
         class EProperties, class GProperties >
class SiconosGraph
{
public:

  /* note : OutEdgeList as multisetS => cannot compile remove_out_edge_if :
     /usr/include/boost/graph/detail/adjacency_list.hpp:440: error: passing 'const ... */
  typedef boost::adjacency_list < boost::listS, boost::listS,
          boost::undirectedS > proxy_graph_t;

  typedef typename
  boost::graph_traits<proxy_graph_t>::edge_descriptor EDescriptor;

  typedef typename
  boost::graph_traits<proxy_graph_t>::vertex_descriptor VDescriptor;


  typedef boost::adjacency_list <
    boost::listS, boost::listS, boost::undirectedS,
    boost::property
    < vertex_siconos_bundle_t, V,
      boost::property <
        boost::vertex_color_t ,
        boost::default_color_type ,
        boost::property < boost::vertex_index_t,
                          size_t,
                          boost::property < vertex_properties_t,
                                            VProperties
                                            > > > > ,
    boost::property
    < edge_siconos_bundle_t, E,
      boost::property <
        boost::edge_color_t ,
        boost::default_color_type ,
        boost::property < boost::edge_index_t,
                          size_t,
                          boost::property < edge_properties_t ,
                                            EProperties
                                            > > > > ,
    boost::property < graph_properties_t, GProperties > >
    graph_t;

  typedef V vertex_t;

  typedef E edge_t;

  typedef typename
  boost::graph_traits<graph_t>::edge_iterator EIterator;

  typedef typename
  boost::graph_traits<graph_t>::vertex_iterator VIterator;

  //  typedef typename
  //  graph_traits<graph_t>::edge_descriptor EDescriptor;

  //  typedef typename
  //  graph_traits<graph_t>::vertex_descriptor VDescriptor;

  typedef typename
  boost::graph_traits<graph_t>::out_edge_iterator OEIterator;

  typedef typename
  boost::graph_traits<graph_t>::adjacency_iterator AVIterator;

  typedef typename
  boost::property_map<graph_t, edge_siconos_bundle_t >::type EBundleAccess;

  typedef typename
  boost::property_map<graph_t, vertex_siconos_bundle_t >::type VBundleAccess;

  typedef typename
  boost::property_map<graph_t, boost::edge_color_t >::type EColorAccess;

  typedef typename
  boost::property_map<graph_t, boost::vertex_color_t >::type VColorAccess;

  typedef typename
  boost::property_map<graph_t, boost::edge_index_t >::type EIndexAccess;

  typedef typename
  boost::property_map<graph_t, boost::vertex_index_t >::type VIndexAccess;

  typedef typename
  boost::property_map<graph_t, edge_properties_t >::type EPropertiesAccess;

  typedef typename
  boost::property_map<graph_t, vertex_properties_t >::type VPropertiesAccess;

  // Need some Base, otherwise, we have a compile error with boost >= 1.51
  //  typedef typename
  //  boost::property_map<graph_t, graph_properties_t >::type GraphPropertiesAccess;

#if defined(SICONOS_STD_UNORDERED_MAP) && !defined(SICONOS_USE_MAP_FOR_HASH)
  typedef typename std::unordered_map<V, VDescriptor> VMap;
#else
  typedef typename std::map<V, VDescriptor> VMap;
#endif


  int _stamp;
  VMap vertex_descriptor;

protected:
  /** serialization hooks
  */
  typedef void serializable;
  template<typename Archive>
  friend void siconos_io(Archive&, SiconosGraph < V, E, VProperties, EProperties,
                         GProperties > &,
                         const unsigned int);
  friend class boost::serialization::access;

  graph_t g;

private:

  SiconosGraph(const SiconosGraph&);

public:

  /** default constructor
   */
  SiconosGraph() : _stamp(0)
  {
  };

  ~SiconosGraph()
  {
    g.clear();
  };

  const graph_t& storage() const
  {
    return g;
  }


  std::pair<EDescriptor, bool>
  edge(VDescriptor u, VDescriptor v) const
  {
    return boost::edge(u, v, g);
  }

  bool edge_exists(const VDescriptor& vd1, const VDescriptor& vd2) const
  {
    bool ret = false;
    EDescriptor tmped;
    std11::tie(tmped, ret) = edge(vd1, vd2);

#ifndef NDEBUG
    bool check_ret = false;
    AVIterator avi, aviend;
    for (std11::tie(avi, aviend) = adjacent_vertices(vd1);
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

  /* parallel edges : edge_range needs multisetS as
     OutEdgesList and with multisetS remove_out_edges_if cannot
     compile as with listS.
     This is only needed for AdjointGraph where only 2 edges may be in
     common which correspond to the source and target in primal graph
   */
  std::pair<EDescriptor, EDescriptor>
  edges(VDescriptor u, VDescriptor v) const
  {
    //    BOOST_STATIC_ASSERT((GProperties::is_adjoint_graph));

    OEIterator oei, oeiend;
    bool ifirst = false;
    bool isecond = false;
    EDescriptor first, second;
    for (std11::tie(oei, oeiend) = out_edges(u); oei != oeiend; ++oei)
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
               const E& e_bundle) const
  {
    bool found = false;
    OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = out_edges(vd1);
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

  bool adjacent_vertex_exists(const VDescriptor& vd) const
  {
    bool ret = false;
    VIterator vi, viend;
    for (std11::tie(vi, viend) = vertices(); vi != viend; ++vi)
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
    return boost::num_vertices(g);
  };

  size_t vertices_number() const
  {
    return boost::num_vertices(g);
  };

  size_t edges_number() const
  {
    return boost::num_edges(g);
  };

  inline V& bundle(const VDescriptor& vd)
  {
    return boost::get(vertex_siconos_bundle, g)[vd];
  };

  inline const V& bundle(const VDescriptor& vd) const
  {
    return boost::get(vertex_siconos_bundle, g)[vd];
  };

  inline E& bundle(const EDescriptor& ed)
  {
    return boost::get(edge_siconos_bundle, g)[ed];
  };
  
  inline const E& bundle(const EDescriptor& ed) const
  {
    return boost::get(edge_siconos_bundle, g)[ed];
  };

  inline boost::default_color_type& color(const VDescriptor& vd)
  {
    return boost::get(boost::vertex_color, g)[vd];
  };

  inline const boost::default_color_type& color(const VDescriptor& vd) const
  {
    return boost::get(boost::vertex_color, g)[vd];
  };

  inline boost::default_color_type& color(const EDescriptor& ed)
  {
    return boost::get(boost::edge_color, g)[ed];
  };

  inline const boost::default_color_type& color(const EDescriptor& ed) const
  {
    return boost::get(boost::edge_color, g)[ed];
  };

  inline GProperties& properties()
  {
    return boost::get_property(g, graph_properties);
  };

  inline size_t& index(const VDescriptor& vd)
  {
    return boost::get(boost::vertex_index, g)[vd];
  };

  inline const size_t& index(const VDescriptor& vd) const
  {
    return boost::get(boost::vertex_index, g)[vd];
  };

  inline size_t& index(const EDescriptor& ed)
  {
    return boost::get(boost::edge_index, g)[ed];
  };

  inline const size_t& index(const EDescriptor& ed) const
  {
    return boost::get(boost::edge_index, g)[ed];
  };

  inline VProperties& properties(const VDescriptor& vd)
  {
    return boost::get(vertex_properties, g)[vd];
  };

  inline EProperties& properties(const EDescriptor& ed)
  {
    return boost::get(edge_properties, g)[ed];
  };

  //  inline VDescriptor& descriptor0(const VDescriptor& vd)
  //  {
  //    return get(vertex_descriptor0, g)[vd];
  //  }

  //  inline EDescriptor& descriptor0(const EDescriptor& ed)
  //  {
  //    return get(edge_descriptor0, g)[ed];
  //  }

  inline bool is_vertex(const V& vertex) const
  {
    return (vertex_descriptor.find(vertex) != vertex_descriptor.end());
  }

  inline const VDescriptor& descriptor(const V& vertex) const
  {
    assert(size() == vertex_descriptor.size());
    assert(vertex_descriptor.find(vertex) != vertex_descriptor.end());
    return (*vertex_descriptor.find(vertex)).second;
  }

  inline std::pair<VIterator, VIterator> vertices() const
  {
    return boost::vertices(g);
  };

  inline VIterator begin() const
  {
    VIterator vi, viend;
    std11::tie(vi, viend) = vertices();
    return vi;
  }

  inline VIterator end() const
  {
    VIterator vi, viend;
    std11::tie(vi, viend) = vertices();
    return viend;
  }

  inline std::pair<AVIterator, AVIterator> adjacent_vertices(const VDescriptor& vd) const
  {
    return boost::adjacent_vertices(vd, g);
  };

  inline std::pair<EIterator, EIterator> edges() const
  {
    return boost::edges(g);
  };

  inline std::pair<OEIterator, OEIterator> out_edges(const VDescriptor& vd) const
  {
    return boost::out_edges(vd, g);
  };

  inline VDescriptor target(const EDescriptor& ed) const
  {
    return boost::target(ed, g);
  };

  inline VDescriptor source(const EDescriptor& ed) const
  {
    return boost::source(ed, g);
  };


  VDescriptor add_vertex(const V& vertex_bundle)
  {
    assert(vertex_descriptor.size() == size()) ;

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

      index(new_vertex_descriptor) = std::numeric_limits<size_t>::max() ;
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
    //    descriptor0(descr) = og.descriptor(vertex_bundle);

    assert(bundle(descr) == vertex_bundle);

    // edges copy as in boost::subgraph
    typename G::OEIterator ogoei, ogoeiend;
    for (std11::tie(ogoei, ogoeiend) =
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

        properties(edescr) = og.properties(*ogoei);
        //        descriptor0(edescr) = *ogoei;

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
    /*   debug */
#ifndef NDEBUG
    assert(adjacent_vertices_ok());
#endif
    boost::clear_vertex(vd, g);
    /*   debug */
#ifndef NDEBUG
    assert(adjacent_vertices_ok());
    assert(!adjacent_vertex_exists(vd));
#endif
    boost::remove_vertex(vd, g);

    assert(vertex_descriptor.size() == (size() + 1));

    vertex_descriptor.erase(vertex_bundle);

    /*  debug */
#ifndef NDEBUG
    assert(adjacent_vertices_ok());
    assert(vertex_descriptor.size() == size());
    assert(!is_vertex(vertex_bundle));
    assert(state_assert());
#endif
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

    std11::tie(new_edge, inserted) = boost::add_edge(vd1, vd2, g);

    // During a gdb session, I saw that inserted is always going to be true ...
    // This check is therefore unnecessary.
    assert(inserted);

    index(new_edge) = std::numeric_limits<size_t>::max();

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

#if defined(SICONOS_STD_UNORDERED_MAP) && !defined(SICONOS_USE_MAP_FOR_HASH)
      std::unordered_map<E, EDescriptor> Edone;
#else
      std::map<E, EDescriptor> Edone;
#endif


      OEIterator ied, iedend;
      for (std11::tie(ied, iedend) = out_edges(vdx);
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
	    //            assert(ag.bundle(ag.add_edge(new_ve, ag.descriptor(bundle(*ied)), bundle(vdx))) == bundle(vdx));
          }
        }

      }
    }
    assert(ag.size() == edges_number());
    return std::pair<EDescriptor, typename AdjointG::VDescriptor>(new_ed,
           new_ve);
  }

  void remove_edge(const EDescriptor& ed)
  {
    assert(adjacent_vertex_exists(target(ed)));
    assert(adjacent_vertex_exists(source(ed)));

    boost::remove_edge(ed, g);
    /* debug */
#ifndef NDEBUG
    assert(state_assert());
#endif
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
    /* debug */
#ifndef NDEBUG
    assert(state_assert());
#endif
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

    BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<graph_t>));
    BOOST_CONCEPT_ASSERT((boost::MutableGraphConcept<graph_t>));

    boost::remove_out_edge_if(vd, pred, g);
    /* workaround on multisetS (tested on Disks : ok)
       multiset allows for member removal without invalidating iterators

        OEIterator oei,oeiend;
        for(std11::tie(oei,oeiend)=out_edges(vd); oei!=oeiend; ++oei)
        {
          if (pred(*oei))
          {
            remove_edge(*oei);
          }
        }
    */
    /*  debug */
#ifndef NDEBUG
    assert(state_assert());
#endif
  }


  /** Remove all the in-edges of vertex u for which the predicate p
   * returns true. This expression is only required when the graph
   * also models IncidenceGraph.
   */
  template<class Predicate>
  void remove_in_edge_if(const VDescriptor& vd,
                          const Predicate& pred)
  //                      Predicate pred)
  {

    BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<graph_t>));
    BOOST_CONCEPT_ASSERT((boost::MutableGraphConcept<graph_t>));

    boost::remove_in_edge_if(vd, pred, g);
    /*  debug */
#ifndef NDEBUG
    assert(state_assert());
#endif
  }

  /** Remove all the in-edges of vertex u for which the predicate p
   * returns true. This expression is only required when the graph
   * also models IncidenceGraph.
   */
  template<class Predicate>
  void remove_edge_if(const VDescriptor& vd,
                      const Predicate& pred)
  //                  Predicate pred)
  {

    BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<graph_t>));
    BOOST_CONCEPT_ASSERT((boost::MutableGraphConcept<graph_t>));

    boost::remove_edge_if(pred, g);
    /*  debug */
#ifndef NDEBUG
    assert(state_assert());
#endif
  }


  int stamp() const
  {
    return _stamp;
  }

  void update_vertices_indices()
  {
    VIterator vi, viend;
    size_t i;
    for (std11::tie(vi, viend) = boost::vertices(g), i = 0;
         vi != viend; ++vi, ++i)
    {
      index(*vi) = i;
    }
    _stamp++;
  };

  void update_edges_indices()
  {
    EIterator ei, eiend;
    size_t i;
    for (std11::tie(ei, eiend) = boost::edges(g), i = 0;
         ei != eiend; ++ei, ++i)
    {
      index(*ei) = i;
    }
    _stamp++;
  };
  
  void clear()
  {
    g.clear();
    vertex_descriptor.clear();
  };

  VMap vertex_descriptor_map() const
  {
    return vertex_descriptor;
  };

  void display() const
  {
    std::cout << "vertices number :" << vertices_number() << std::endl;
        
    std::cout << "edges number :" << edges_number() << std::endl;
    VIterator vi, viend;
    for (std11::tie(vi, viend) = vertices();
         vi != viend; ++vi)
    {
      std::cout << "vertex :"
                << *vi
                << ", bundle :"
                << bundle(*vi)
                << ", index : "
                << index(*vi)
                << ", color : "
                << color(*vi);
      OEIterator oei, oeiend, next;
      for (std11::tie(oei, oeiend) = out_edges(*vi);
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
#ifndef SWIG
#ifndef NDEBUG
  bool state_assert() const
  {
    VIterator vi, viend;
    for (std11::tie(vi, viend) = vertices(); vi != viend; ++vi)
    {
      assert(is_vertex(bundle(*vi)));
      assert(bundle(descriptor(bundle(*vi))) == bundle(*vi));

      OEIterator ei, eiend;
      for (std11::tie(ei, eiend) = out_edges(*vi);
           ei != eiend; ++ei)
      {
        assert(is_vertex(bundle(target(*ei))));
        assert(source(*ei) == *vi);
      }
      AVIterator avi, aviend;
      for (std11::tie(avi, aviend) = adjacent_vertices(*vi);
           avi != aviend; ++avi)
      {
        assert(is_vertex(bundle(*avi)));
        assert(bundle(descriptor(bundle(*avi))) == bundle(*avi));
      }
    }
    return true;

  }

  bool adjacent_vertices_ok() const
  {

    VIterator vi, viend;
    for (std11::tie(vi, viend) = vertices(); vi != viend; ++vi)
    {
      assert(is_vertex(bundle(*vi)));
      assert(bundle(descriptor(bundle(*vi))) == bundle(*vi));

      AVIterator avi, aviend;
      for (std11::tie(avi, aviend) = adjacent_vertices(*vi);
           avi != aviend; ++avi)
      {
        assert(is_vertex(bundle(*avi)));
        assert(bundle(descriptor(bundle(*avi))) == bundle(*avi));
      }
    }
    return true;
  }
#endif
#endif


};


#endif


