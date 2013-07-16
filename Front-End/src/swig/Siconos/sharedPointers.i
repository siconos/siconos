%{
#include <SiconosPointers.hpp>
#include <SiconosFwd.hpp>
%}
#define SWIG_SHARED_PTR_NAMESPACE std11
%include boost_shared_ptr.i

#if (__cplusplus >= 201103L) && !defined(USE_BOOST_FOR_CXX11)
#define STD11 std
#else
#define STD11 boost
%import "boost/version.hpp"
//  boost >= 1.53
// this sucks and will likely not work in C++11, but it is difficult to
// deal with this properly
#if (BOOST_VERSION >= 105300)
#define BOOST_NOEXCEPT
#endif
// boost >= 1.40
#if (BOOST_VERSION >= 104000)
%ignore std11::enable_shared_from_this::operator=;
%import "boost/smart_ptr/enable_shared_from_this.hpp"
#else
%import "boost/enable_shared_from_this.hpp"
#endif
#endif

%rename("$ignore", regexmatch$name="^createSPtr.*") "";

%{
  // when we call FPyArray_SimpleNewFromData with a $1->getArray() we
  // lost the shared pointer count, so we can be in the case where the
  // memory pointed by a shared ptr is erased after the call
  // FPyArray_SimpleNewFromData (=>segfault...)

  // here we keep another shared pointer on the original shared ptr
  // (i.e. we do an incref) associated with the memory from
  // FPyArray_SimpleNewFromData

  // we need to register a PyCObject (deprecated for python >= 2.7 see
  // PyCapsule) in order to call the destruction function
  // sharedPyArrayDelete


  struct SharedPointerKeeper
  {
    // to keep a pointer on shared_ptr{Siconos,Simple}{Vector,Matrix}
    std11::shared_ptr<void> ref;

    SharedPointerKeeper(std11::shared_ptr<void> v) : ref(v) 
    {
      DEBUG_PRINTF("SharedPointerKeeper : use_count %ld\n",v.use_count());
    };

    ~SharedPointerKeeper()
    {
      DEBUG_PRINT("~SharedPointerKeeper()\n");
      //    ref.reset(); // destructor called
    }

  };
  
  /* the PyCObject deleter 
     example: 
     SharedPointerKeeper* savedSharePtr = 
       new SharedPointerKeeper(std11::static_pointer_cast<void>(mysharedptr));
     PyCObject_FromVoidPtr((void*) savedSharedPtr, &sharedPointerKeeperDelete);
  */

  /* note PyCObject is deprecated for Python >= 2.7 ... */
  static  void sharedPointerKeeperDelete(void * o)
  {
    DEBUG_PRINT("sharedPointerKeeperDelete\n");

    delete static_cast<SharedPointerKeeper *>(o);
    return;
  };

%}


%include SiconosPointers.hpp
%include SiconosFwd.hpp
