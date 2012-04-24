// Copyright (C) Vladimir Prus 2003.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/graph/vector_property_map.html for
// documentation.
//

#ifndef SICONOS_PROPERTIES_HPP
#define SICONOS_PROPERTIES_HPP

#include <boost/property_map/property_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/static_assert.hpp>
#include <vector>

namespace Siconos
{


/* a templated way to access edge or vertex with the help of
 * boost::mpl::if_ */

template<typename G>
struct VertexAccess
{
  typedef VertexAccess type;
  typedef typename G::VDescriptor descriptor;
  typedef typename G::OVIndexAccess OldIndexMap;

  std::vector<descriptor>& indexModified(G& g)
  {
    return g.vertex_index_modified();
  }

  size_t size(G& g)
  {
    return g.vertices_number();
  }

};

template<typename G>
struct EdgeAccess
{
  typedef EdgeAccess type;
  typedef typename G::EDescriptor descriptor;
  typedef typename G::OEIndexAccess OldIndexMap;

  static std::vector<descriptor>& indexModified(G& g)
  {
    return g.edge_index_modified();
  }

  static size_t size(const G& g)
  {
    return g.edges_number();
  }
};


template<typename G, typename IndexMap>
struct VertexOrEdge
{
  typedef typename boost::mpl::if_ < boost::is_same<typename G::VIndexAccess, IndexMap>,
          VertexAccess<G> ,
          EdgeAccess<G>  >::type Access;
};


/* the properties structure, from boost::vector_property_map */

template<typename T, typename G, typename IndexMap>
class Properties
{

private:

  boost::shared_ptr<G> _g;
  boost::shared_ptr< std::vector<T> > _store;
  int _stamp;

public:
  typedef typename VertexOrEdge<G, IndexMap>::Access Access;

  Access access;

  typedef typename Access::OldIndexMap OldIndexMap;
  typedef typename boost::property_traits<IndexMap>::key_type  key_type;
  typedef T value_type;
  typedef typename std::iterator_traits <
  typename std::vector<T>::iterator >::reference reference;
  typedef boost::lvalue_property_map_tag category;

  Properties(boost::shared_ptr<G> g)
    : _g(g), _store(new std::vector<T>()), _stamp(0)
  {}

  reference operator[](const key_type& v)
  {
    if (_g->stamp() != _stamp)
    {
      _stamp = _g->stamp();

      int k = access.size(*_g);
      typename std::vector<T> new_store(k);

      k--;
      for (typename std::vector<typename Access::descriptor>::iterator vi = access.indexModified(*_g).begin();
           vi != access.indexModified(*_g).end(); ++vi)
      {
        typename boost::property_traits<OldIndexMap>::value_type oi = _g->old_index(*vi);
        new_store[k--] = (*_store)[oi];
      }
      for (typename std::vector<typename Access::descriptor>::iterator vi = access.indexModified(*_g).begin();
           vi != access.indexModified(*_g).end(); ++vi)
      {
        typename boost::property_traits<IndexMap>::value_type i = _g->index(*vi);
        (*_store)[i] = new_store.back();
        new_store.pop_back();
      };
    };

    typename boost::property_traits<IndexMap>::value_type i = _g->index(v);
    if (static_cast<unsigned>(i) >= _store->size())
    {
      // _store->resize(i+1, T()); faster with large size
      _store->resize(access.size(*_g), T());
    }

    return (*_store)[i];
  };

};

template<typename T, typename G>
class VertexProperties : public Properties<T, G, typename G::VIndexAccess>
{
public:
  VertexProperties(boost::shared_ptr<G> g) : Properties<T, G, typename G::VIndexAccess>(g)
  {};
};

template<typename T, typename G>
class EdgeProperties : public Properties<T, G, typename G::EIndexAccess>
{
public:
  EdgeProperties(boost::shared_ptr<G> g) : Properties<T, G, typename G::EIndexAccess>(g)
  {};
};


template<typename T, typename G>
VertexProperties<T, G> vertexProperties(boost::shared_ptr<G> g)
{
  return VertexProperties<T, G>(g);
}

template<typename T, typename G>
EdgeProperties<T, G> edgeProperties(boost::shared_ptr<G> g)
{
  return EdgeProperties<T, G>(g);
}

}

/* convenience macro for properties declarations */

#include <boost/preprocessor/seq/seq.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/cat.hpp>

#define I_DECLARE_MEMBERS(r,gt,p) \
  Siconos:: BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(3,0,p),Properties)< BOOST_PP_TUPLE_ELEM(3,1,p), gt> BOOST_PP_TUPLE_ELEM(3,2,p);

#define I_CONS_MEMBERS(r,gt,p) \
  BOOST_PP_TUPLE_ELEM(3,2,p) ( Siconos:: BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(3,0,p),Properties)< BOOST_PP_TUPLE_ELEM(3,1,p), gt>(g) ),


#define INSTALL_PROPERTIES(GraphType, PROPERTIES)                       \
  struct BOOST_PP_CAT(GraphType,Properties) : public GraphProperties    \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(I_DECLARE_MEMBERS, GraphType, PROPERTIES);    \
                                                                        \
    bool dummy;                                                         \
                                                                        \
    GraphType##Properties(SP:: GraphType g) :                           \
      BOOST_PP_SEQ_FOR_EACH(I_CONS_MEMBERS, GraphType, PROPERTIES) dummy(false) {}; \
  }


#endif
