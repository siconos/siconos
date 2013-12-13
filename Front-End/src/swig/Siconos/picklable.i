
#ifndef PICKLABLE_i
#define PICKLABLE_i

#ifdef WITH_IO
%{
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <sstream>
%}

/* allow python serialization from SiconosIO serializers */
%define %make_picklable(CLASS)
%extend CLASS {
  std::string __getstate__()
  {
    std::stringstream ss;
    boost::archive::binary_oarchive ar(ss);
    ar << *($self);
    return ss.str();
  }
  
  void __setstate_internal(std::string const& from_str)
  {
    std::stringstream ss(from_str);
    boost::archive::binary_iarchive ar(ss);
    ar >> *($self);
  }
  
  %pythoncode %{
    def __setstate__(self, from_str):
        self.__init__()
        self.__setstate_internal(from_str)
      %}

  std::string xml_export()
  {
    std::stringstream ss;
    boost::archive::xml_oarchive ar(ss);
    ar << BOOST_SERIALIZATION_NVP(*($self));
    return ss.str();
  }

  void xml_import(std::string from_str)
  {
    std::stringstream ss(from_str);
    boost::archive::xml_iarchive ar(ss);
    ar >> BOOST_SERIALIZATION_NVP(*($self));
  }
  
  std::string text_export()
  {
    std::stringstream ss;
    boost::archive::text_oarchive ar(ss);
    ar << BOOST_SERIALIZATION_NVP(*($self));
    return ss.str();
  }
  
  void text_import(std::string from_str)
  {
    std::stringstream ss(from_str);
    boost::archive::text_iarchive ar(ss);
    ar >> BOOST_SERIALIZATION_NVP(*($self));
  }
  
 }
%enddef
#endif

#ifndef WITH_IO

%define %make_picklable(CLASS)
%enddef

#endif


#endif
