
#ifdef WITH_SERIALIZATION
%{
#include <SiconosFull.hpp>

// Work-around for issue reading inf/nan double values
#include <boost/archive/basic_text_iprimitive.hpp>
namespace boost { namespace archive {
template<> template<>
void basic_text_iprimitive<std::istream>::load<double>( double& t )
{
  char s[32]={}, *b, *a = b = s;
  int i=0;
  for (int i=0; is.peek()!='<' && i < 31; i++)
    s[i] = is.get();
  errno = 0;
  t = std::strtod(s, &b);
  if (errno == ERANGE || errno == EINVAL || a==b)
    boost::serialization::throw_exception(
      archive_exception(archive_exception::input_stream_error));
}
} } // close namespaces
%}
#endif
