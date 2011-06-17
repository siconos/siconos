#ifndef MEMBERS_HPP
#define MEMBERS_HPP

#include <boost/preprocessor/seq/seq.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

#include <boost/mpl/eval_if.hpp>


#define INTERNAL_SICONOS_SERIALIZATION_NVP(class,member)                \
  boost::serialization::make_nvp(BOOST_PP_STRINGIZE(member), class.member)

/* serialization is not splitted */
#define INTERNAL_SICONOS_IO_SERIALIZE(r,c,m) \
  ar & INTERNAL_SICONOS_SERIALIZATION_NVP(c,m);

#define INTERNAL_SICONOS_IO_SERIALIZE_BASE(r,c,b)                       \
  ar & boost::serialization::make_nvp(                                  \
    BOOST_PP_STRINGIZE(b),                                              \
    boost::serialization::base_object<b>(c) );                          \
 
#define INTERNAL_SICONOS_BOOST_REGISTER(CLASS)                          \
  namespace boost { namespace serialization {                           \
      struct Load##CLASS                                                \
      {                                                                 \
        typedef Load##CLASS type;                                       \
        template<typename Archive>                                      \
          void function(Archive& ar, CLASS& c, const unsigned int version) \
        { load(ar,c,version); };                                        \
      };                                                                \
      struct Save##CLASS                                                \
      {                                                                 \
        typedef Save##CLASS type;                                       \
        template<typename Archive>                                      \
          void function(Archive& ar, CLASS& c, const unsigned int version) \
        { save(ar,c,version); };                                        \
      };                                                                \
      template<class Archive>                                           \
      void serialize(Archive & ar, CLASS & c, const unsigned int version) \
      {                                                                 \
        typedef typename boost::mpl::eval_if<typename Archive::is_saving, \
          Load##CLASS, Save##CLASS>::type Serializer;                   \
        static Serializer serializer;                                   \
        serializer.function(ar,c,version);};                            \
    }}

/** base class members registration
    \param class name
    \param members sequence (as a boost sequence (member1)(member2) ...)
*/
#define SICONOS_IO_REGISTER(CLASS,MEMBERS)                              \
  template<class Archive>                                               \
  void save(Archive & ar, CLASS & c, const unsigned int version)        \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, c, MEMBERS);   \
  };                                                                    \
  template<class Archive>                                               \
  void load(Archive & ar, CLASS & c, const unsigned int version)        \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, c, MEMBERS);   \
  };                                                                    \
  INTERNAL_SICONOS_BOOST_REGISTER(CLASS)



/** derived class members registration
    \param derived class name
    \param base class name
    \param members sequence (as a boost sequence (member1)(member2) ...)
*/
#define SICONOS_IO_REGISTER_WITH_BASE(CLASS,BASE,MEMBERS)               \
  template<class Archive>                                               \
  void save(Archive & ar, CLASS & c, const unsigned int version)        \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, c, MEMBERS);   \
    ar << boost::serialization::make_nvp(                               \
      BOOST_PP_STRINGIZE(CLASS),                                        \
      boost::serialization::base_object<BASE>(c) );                     \
  };                                                                    \
  template<class Archive>                                               \
  void load(Archive & ar, CLASS & c, const unsigned int version)        \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, c, MEMBERS);   \
    ar >> boost::serialization::make_nvp(                               \
      BOOST_PP_STRINGIZE(BASE),                                         \
      boost::serialization::base_object<BASE>(c) );                     \
  };                                                                    \
  INTERNAL_SICONOS_BOOST_REGISTER(CLASS)


/** derived class with multiple inheritance registration
    \param derived class name
    \param base class name
    \param members sequence (as a boost sequence (member1)(member2) ...)
*/
#define SICONOS_IO_REGISTER_WITH_BASES(CLASS,BASES,MEMBERS)             \
  template<class Archive>                                               \
  void save(Archive & ar, CLASS & c, const unsigned int version)        \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, c, MEMBERS);   \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE_BASE, c, BASES); \
  };                                                                    \
  template<class Archive>                                               \
  void load(Archive & ar, CLASS & c, const unsigned int version)        \
  {                                                                     \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE, c, MEMBERS);   \
    BOOST_PP_SEQ_FOR_EACH(INTERNAL_SICONOS_IO_SERIALIZE_BASE, c, BASES); \
  };                                                                    \
  INTERNAL_SICONOS_BOOST_REGISTER(CLASS)



#endif
