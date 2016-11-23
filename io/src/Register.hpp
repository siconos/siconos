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

#ifndef Register_hpp
#define Register_hpp

#include "SiconosConfig.h"

#include <boost/preprocessor/seq/seq.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

#if defined(SICONOS_STD_SHARED_PTR) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
#include <boost/serialization/ser_shared_ptr.hpp>
#else
#include <boost/serialization/shared_ptr.hpp>
#endif

#include <boost/serialization/weak_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/hash_set.hpp>
#include <boost/serialization/deque.hpp>
#include "boost/serialization/unordered_set.hpp"

#include <boost/serialization/list.hpp>

#include <boost/graph/adj_list_serialize.hpp>

#include <boost/serialization/export.hpp>

#include <boost/typeof/typeof.hpp>

#include <debug.h>

// This should be cleaner --xhub
#if (BOOST_VERSION >= 106100)
#define boost_ser_array boost::serialization::array_wrapper
#else
#define boost_ser_array boost::serialization::array
#endif

/** register class serialization in boost namespace
    \param a class name
 */
#define REGISTER_BOOST_SERIALIZATION(C)                                 \
  namespace boost { namespace serialization {                           \
      template <class Archive>                                          \
      void serialize(Archive& ar, C& v, unsigned int version)           \
      {                                                                 \
        siconos_io(ar, v, version);                                     \
      }                                                                 \
    }                                                                   \
  }                                                                     \

/** internal macro */
#define INTERNAL_SICONOS_SERIALIZATION_NVP(object,member)               \
  ::boost::serialization::make_nvp(BOOST_PP_STRINGIZE(member), object.member)

/** internal macro */
#define INTERNAL_SICONOS_IO_SERIALIZE(r,o,m) \
  ar & INTERNAL_SICONOS_SERIALIZATION_NVP(o,m);

/** internal macro */
#define INTERNAL_SICONOS_IO_SERIALIZE_BASE(r,o,b)                       \
  ar & ::boost::serialization::make_nvp(                                \
    BOOST_PP_STRINGIZE(b),                                              \
    ::boost::serialization::base_object<b>(o) );                        \
 
/** base class members registration
 *  \param class name
 *  \param members sequence (as a boost preprocessor sequence
 *   (member1)(member2)x...)
 */
#define SICONOS_IO_REGISTER(CLASS,MEMBERS)                              \
  template<class Archive>                                               \
  void siconos_io(Archive & ar, CLASS & o, const unsigned int version) \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, o, MEMBERS);   \
  };                                                                    \
  namespace boost { namespace serialization {                           \
      template<class Archive>                                           \
      void serialize(Archive & ar, CLASS & o, const unsigned int version) \
      {                                                                 \
        DEBUG_PRINTF("serialize %s\n", BOOST_PP_STRINGIZE(CLASS));      \
        siconos_io(ar,o,version);                                       \
      };                                                                \
    }}

/** derived class members registration
    \param derived class name
    \param base class name
    \param members sequence (as a boost sequence (member1)(member2) ...)
*/
#define SICONOS_IO_REGISTER_WITH_BASE(CLASS,BASE,MEMBERS)               \
  template<class Archive>                                               \
  void siconos_io(Archive & ar, CLASS & o, const unsigned int version)  \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, o, MEMBERS);   \
    ar &  ::boost::serialization::make_nvp(                             \
      BOOST_PP_STRINGIZE(BASE),                                         \
      ::boost::serialization::base_object<BASE>(o) );                   \
  };                                                                    \
  namespace boost { namespace serialization {                           \
      template<class Archive>                                           \
      void serialize(Archive & ar, CLASS & o, const unsigned int version) \
      {                                                                 \
        DEBUG_PRINTF("serialize %s\n", BOOST_PP_STRINGIZE(CLASS));      \
        siconos_io(ar,o,version);                                       \
      };                                                                \
    }}

/** derived class with multiple inheritance registration
    \param derived class name
    \param base class name
    \param members sequence (as a boost sequence (member1)(member2) ...)
*/
#define SICONOS_IO_REGISTER_WITH_BASES(CLASS,BASES,MEMBERS)             \
  template<class Archive>                                               \
  void siconos_io(Archive & ar, CLASS & o, const unsigned int version)  \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, o, MEMBERS);   \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE_BASE, o, BASES); \
  };                                                                    \
  namespace boost { namespace serialization {                           \
      template<class Archive>                                           \
      void serialize(Archive & ar, CLASS & o, const unsigned int version) \
      {                                                                 \
        siconos_io(ar,o,version);                                       \
      };                                                                \
    }}

#include <boost/preprocessor/tuple/elem.hpp>

/** internal macro */
#define SERIALIZE_I(r,T,M)                                              \
  BOOST_PP_TUPLE_ELEM(2,0,T) & \
  ::boost::serialization::make_nvp(BOOST_PP_STRINGIZE(M), \
                                   BOOST_PP_TUPLE_ELEM(2,1,T) . M);

/** serialize structure members
 * \param a structure instance
 * \param a boost preprocessor sequence of members
 * \param an archive
 */
#define SERIALIZE(S,MEMBERS, ARCHIVE)                       \
  BOOST_PP_SEQ_FOR_EACH(SERIALIZE_I, (ARCHIVE, S), MEMBERS)

/** serialize C array inside structure
 * \param array dimension
 * \param a struct instance
 * \param array member
 * \param an archive
 */
#define SERIALIZE_C_ARRAY(DIM, STRUCT, ARRAY, ARCHIVE)                   \
  if (Archive::is_loading::value)                                       \
  {                                                                     \
    STRUCT . ARRAY = (BOOST_TYPEOF(STRUCT . ARRAY)) malloc(DIM * sizeof(BOOST_TYPEOF(* (STRUCT . ARRAY)))); \
  };                                                                    \
  {                                                                     \
    boost_ser_array<BOOST_TYPEOF(*(STRUCT . ARRAY))>        \
      wrapper = boost::serialization::make_array(STRUCT . ARRAY,DIM);   \
    ARCHIVE & boost::serialization::make_nvp(BOOST_PP_STRINGIZE(ARRAY),wrapper); \
  }


#endif
