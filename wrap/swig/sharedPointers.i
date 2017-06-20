//%include "SiconosConfig.h"
%{
#include <SiconosPointers.hpp>
#include <SiconosFwd.hpp>
%}

#if defined(SICONOS_STD_SHARED_PTR) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
#define STD11 std
#undef __cplusplus
#define __cplusplus SICONOS_CXXVERSION
%include <std_shared_ptr.i>

// from g++-v4/bits/shared_ptr.h
// not sure if this is needed, but we can't use '#include <memory>'
// since it is in a compiler path
namespace STD11 {
  template<typename _Tp>
    class enable_shared_from_this
  {
  protected:
    constexpr enable_shared_from_this();

    ~enable_shared_from_this();

  public:
    shared_ptr<_Tp>
      shared_from_this();

    shared_ptr<const _Tp>
      shared_from_this() const;
  };
 }
#else
#define SWIG_SHARED_PTR_NAMESPACE std11
%include <boost_shared_ptr.i>
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
// boost >= 1.64
#if (BOOST_VERSION >= 106400)
%import "boost/smart_ptr/detail/sp_noexcept.hpp"
#endif
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

#ifdef SWIGPY_USE_CAPSULE

  static  void sharedPointerKeeperDeleteCap(PyObject * cap)
  {
    DEBUG_PRINT("sharedPointerKeeperDeleteCap\n");
    void* o = (void*) PyCapsule_GetPointer(cap,SWIGPY_CAPSULE_NAME);
    delete static_cast<SharedPointerKeeper *>(o);
    return;
  };

#else
  /* note PyCObject is deprecated for Python >= 2.7 ... */
  static  void sharedPointerKeeperDelete(void * o)
  {
    DEBUG_PRINT("sharedPointerKeeperDelete\n");

    delete static_cast<SharedPointerKeeper *>(o);
    return;
  };

#endif


%}


%include SiconosPointers.hpp
%include SiconosFwd.hpp
