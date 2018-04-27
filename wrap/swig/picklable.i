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

/* the binary_import/export functions below need to use the 'bytes'
 * type in Python 3, otherwise SWIG tries and fails to convert them to
 * unicode, so we define typemaps for a 'bytes' typedef that
 * distinguishes this from text strings */
typedef std::string bytes;
%{
typedef std::string bytes;
%}

%typemap(out) bytes %{
#if PY_VERSION_HEX >= 0x03000000
  $result = PyBytes_FromStringAndSize($1.c_str(), $1.size());
#else
  $result = PyString_FromStringAndSize($1.c_str(), $1.size());
#endif
%}

%typemap(typecheck) bytes %{
#if PY_VERSION_HEX>=0x03000000
  $1 = PyBytes_Check(obj) ? 1 : 0;
#else
  $1 = PyString_Check(obj) ? 1 : 0;
#endif
%}

%typemap(in) (bytes const&) %{
  {
    char *cstr; Py_ssize_t len;
#if PY_VERSION_HEX>=0x03000000
    PyBytes_AsStringAndSize($input, &cstr, &len);
#else
    PyString_AsStringAndSize($input, &cstr, &len);
#endif
    $1 = new std::string(cstr, cstr+len);
  }
%}

%typemap(freearg) (bytes const&) %{
  if ($1) delete $1;
%}

/* allow python serialization from SiconosIO serializers */
%define %make_picklable(CLASS, COMPONENT)
%extend CLASS {
  bytes binary_export()
  {
    std::stringstream ss;
    boost::archive::binary_oarchive ar(ss);
    siconos_io_register_ ## COMPONENT(ar);
    ar << ::boost::serialization::make_nvp(BOOST_PP_STRINGIZE(CLASS),(*($self)));
    return ss.str();
  }

  bytes __getstate__()
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

  void binary_import(bytes const& from_bytes)
  {
    std::stringstream ss(from_bytes);
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
