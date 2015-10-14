#ifndef BOOST_SERIALIZATION_STD_SHARED_PTR_HPP
#define BOOST_SERIALIZATION_STD_SHARED_PTR_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// shared_ptr.hpp: serialization for std shared pointer

// (C) Copyright 2004 Robert Ramey and Martin Ecker
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <cstddef> // NULL

#include <boost/serialization/split_free.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// shared_ptr serialization traits

namespace boost
{
namespace serialization
{

template<class T>
struct version< ::std::shared_ptr< T > >
{
  static const int value = 0;
};

// don't track shared pointers
template<class T>
struct tracking_level< ::std::shared_ptr< T > >
{

  static const int value = ::boost::serialization::track_never;
};

}
}

namespace boost
{
namespace serialization
{

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// serialization for shared_ptr


template<class Archive, class T>
inline void save(
  Archive & ar,
  const std::shared_ptr< T > &t,
  const unsigned int /* file_version */
)
{
  // The most common cause of trapping here would be serializing
  // something like shared_ptr<int>.  This occurs because int
  // is never tracked by default.  Wrap int in a trackable type
  BOOST_STATIC_ASSERT((tracking_level< T >::value != track_never));
  const T * t_ptr = t.get();
  ar << boost::serialization::make_nvp("px", t_ptr);

  // TODO: cuidado com acesso concorrente!
  static boost::shared_ptr<std::map<T *, std::shared_ptr<T>>> mapping;
  if (!mapping)
    mapping.reset(new std::map<T*, std::shared_ptr<T>>());

  ar << boost::serialization::make_nvp("map", mapping);
}

namespace detail
{
template < class T, class EN =
typename std::enable_if<T::has_shared_from_this>::type >
std::shared_ptr<T> get_smart_ptr(T *ptr)
{
  return std::static_pointer_cast<T>(ptr->shared_from_this());
}

template <class T>
std::shared_ptr<T> get_smart_ptr(void *ptr)
{
  return std::shared_ptr<T>((T *)ptr);
}
}

template<class Archive, class T>
inline void load(
  Archive & ar,
  std::shared_ptr< T > &t,
  const unsigned int /*file_version*/
)
{
  // The most common cause of trapping here would be serializing
  // something like shared_ptr<int>.  This occurs because int
  // is never tracked by default.  Wrap int in a trackable type
  BOOST_STATIC_ASSERT((tracking_level< T >::value != track_never));
  T* r = NULL;
  ar >> boost::serialization::make_nvp("px", r);

  boost::shared_ptr<std::map<T *, std::shared_ptr<T>>> mapping;

  ar >> boost::serialization::make_nvp("map", mapping);

  if (r == NULL)
    t.reset();
  else
  {
    auto it = mapping->find(r);
    if (it == mapping->end())
    {
      std::shared_ptr<T> p = detail::get_smart_ptr<T>(r);
      it = mapping->insert( {r, p}).first;
    }

    t = it->second;
  }
}

template<class Archive, class T>
inline void serialize(
  Archive & ar,
  std::shared_ptr< T > &t,
  const unsigned int file_version
)
{
  // correct shared_ptr serialization depends upon object tracking
  // being used.
  BOOST_STATIC_ASSERT(
    boost::serialization::tracking_level< T >::value
    != boost::serialization::track_never
  );
  boost::serialization::split_free(ar, t, file_version);
}

template<class Archive, class T>
inline void save(
  Archive & ar,
  const std::weak_ptr< T > &t,
  const unsigned int /* file_version */
)
{
  const std::shared_ptr< T > sp = t.lock();
  ar << boost::serialization::make_nvp("weak_ptr", sp);
}

template<class Archive, class T>
inline void load(
  Archive & ar,
  std::weak_ptr< T > &t,
  const unsigned int /* file_version */
)
{
  std::shared_ptr< T > sp;
  ar >> boost::serialization::make_nvp("weak_ptr", sp);
  t = sp;
}

template<class Archive, class T>
inline void serialize(
  Archive & ar,
  std::weak_ptr< T > &t,
  const unsigned int file_version
)
{
  boost::serialization::split_free(ar, t, file_version);
}


} // namespace serialization
} // namespace boost

#endif // BOOST_SERIALIZATION_STD_SHARED_PTR_HPP
