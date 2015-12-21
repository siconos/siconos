#ifndef PICKLABLE_i
#define PICKLABLE_i

#ifdef WITH_SERIALIZATION
%{
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <sstream>
%}

/* allow python serialization from SiconosIO serializers */
%define %make_picklable(CLASS, COMPONENT)
%extend CLASS {
  std::string binary_export()
  {
    std::stringstream ss;
    boost::archive::binary_oarchive ar(ss);
    siconos_io_register_ ## COMPONENT(ar);
    ar << ::boost::serialization::make_nvp(BOOST_PP_STRINGIZE(CLASS),(*($self)));
    return ss.str();
  }

   std::string __getstate__()
  {
    return CLASS##_binary_export($self);
  }

  %pythoncode %{
    def __setstate__(self, from_str):
        self.__init__()
        self.binary_import(from_str)
      %}

  std::string xml_export()
  {
    std::stringstream ss;
    boost::archive::xml_oarchive ar(ss);
    siconos_io_register_ ## COMPONENT(ar);
    ar << ::boost::serialization::make_nvp(BOOST_PP_STRINGIZE(CLASS),(*($self)));
    return ss.str();
  }

  void xml_import(std::string const& from_str)
  {
    std::stringstream ss(from_str);
    boost::archive::xml_iarchive ar(ss);
    siconos_io_register_ ## COMPONENT(ar);
    ar >> ::boost::serialization::make_nvp(BOOST_PP_STRINGIZE(CLASS),(*($self)));
  }


  void binary_import(std::string const& from_str)
  {
    std::stringstream ss(from_str);
    boost::archive::binary_iarchive ar(ss);
    siconos_io_register_ ## COMPONENT(ar);
    ar >> ::boost::serialization::make_nvp(BOOST_PP_STRINGIZE(CLASS),(*($self)));
  }

 }
%enddef
#endif

#ifndef WITH_SERIALIZATION

%define %make_picklable(CLASS, COMPONENT)
%enddef

#endif


#endif
