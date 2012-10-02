// -*- C++ -*-
// Siconos-Front-End , Copyright INRIA 2005-2012.
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.	
// Siconos is a free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// Siconos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Siconos; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr 
//	

// SWIG interface for Siconos Kernel
%module(directors="1", allprotected="1") Kernel

%feature("director:except") {
  if ($error != NULL) {
    throw Swig::DirectorMethodException();
  }
 }


 // signatures
%feature("autodoc", 1);

// generated docstrings from doxygen xml output
%include Kernel-docstrings.i

%{
#define SWIG_FILE_WITH_INIT
#include <sstream>
#if defined(Py_COMPLEXOBJECT_H)
#undef c_sum
#undef c_diff
#undef c_neg
#undef c_prod
#undef c_quot
#undef c_pow
#undef c_abs
#endif
#include <SiconosKernel.hpp>
#include <SiconosVisitor.hpp>
#include "SiconosPointers.hpp"
#include "Disk.hpp"
#include "Circle.hpp"
#include "DiskDiskR.hpp"
#include "DiskPlanR.hpp"
#include "DiskMovingPlanR.hpp"
#include "CircleCircleR.hpp"
#include "ExternalBody.hpp"
#include "SpaceFilter.hpp"
#include "SiconosBodies.hpp"
#include "addons.hpp"
#include "KneeJointR.hpp"
#include "PivotJointR.hpp"
#include "PrismaticJointR.hpp"

#include <boost/type_traits/is_polymorphic.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/typeof/typeof.hpp>

#include <assert.h>

//#define DEBUG_MESSAGES 1
#include <debug.h>

#include "FrontEndConfig.h"
#ifdef HAVE_SICONOS_IO
#include <SiconosRestart.hpp>
#endif
%} 

#ifdef WITH_BULLET
%include "KernelBullet.i"
#endif

// common declarations with Numerics

%include Common.i

// mandatory !
%rename (lambda_) lambda;

// shared ptr management
#define SWIG_SHARED_PTR_NAMESPACE std11
%include "boost_shared_ptr.i"

// numpy macros
%include numpy.i 	

%init %{
  import_array();
%}

// handle standard exceptions
%{
 static void handle_exception(void) {
    try {
      throw;
    } 
    catch (const std::invalid_argument& e)
    {
      PyErr_SetString(PyExc_ValueError, e.what());
    }
    catch (const std::out_of_range& e)
    {
      PyErr_SetString(PyExc_IndexError, e.what());
    }
    catch (const SiconosException& e)
    {
      PyErr_SetString(PyExc_Exception, e.report().c_str());
    }
    catch (const Swig::DirectorException& e)
    {
      PyErr_SetString(PyExc_ValueError, e.getMessage());
    }
  }
%} 


%include "exception.i"
%exception
{
  try
  {
    $action;
  }
  catch (...) {
    if (!PyErr_Occurred()) {
      handle_exception();
    }
    SWIG_fail;
  } 
}

// handle stl data types
%include "stl.i"

// 1. Vector and Matrix <=> numpy array (dense only)

%include "KernelTypes.i"




// python int sequence => std::vector<unsigned int>
%{
  static int sequenceToUnsignedIntVector(
    PyObject *input, 
    std11::shared_ptr<std::vector<unsigned int> > ptr) 
  {
    int i;
    if (!PySequence_Check(input)) {
      PyErr_SetString(PyExc_TypeError,"Expecting a sequence");
      return 0;
    }
    
    assert(ptr);
    
    for (i =0; i <  PyObject_Length(input); i++) 
    {
      PyObject *o = PySequence_GetItem(input,i);
      if (!PyInt_Check(o)) {
        Py_XDECREF(o);
        PyErr_SetString(PyExc_ValueError,"Expecting a sequence of ints");
        return 0;
      }
      
      if (PyInt_AsLong(o) == -1 && PyErr_Occurred())
        return 0;
      
      ptr->push_back(static_cast<unsigned int>(PyInt_AsLong(o)));
      
      Py_DECREF(o);
    }
    return 1;
  }
%}

// 2. try to hide SP::Type on python side


// boost namespace (can be fixed with a correct import)
namespace boost 
{
  namespace numeric
  {
    namespace ublas
    {
    }
    namespace bindings
    {
      namespace atlas
      {
      }
    }
  }
};

// a std::size_t definition (otherwise swig complains about it)
namespace std 
{
  typedef size_t size_t;
}

#if (__cplusplus >= 201103L) && !defined(USE_BOOST_FOR_CXX11)
#define STD11 std
#else
#define STD11 boost
// boost >= 1.40
%import "boost/version.hpp"
#if (BOOST_VERSION >= 104000)
%ignore std11::enable_shared_from_this::operator=;
%import "boost/smart_ptr/enable_shared_from_this.hpp"
#else
%import "boost/enable_shared_from_this.hpp"
#endif
#endif

// SimpleMatrix operators 
%rename  (__add__) operator+;
%rename  (__less__) operator-;
%rename  (__mul__) operator*;
%rename  (__div__) operator/;
%rename  (__iadd__) operator+=;
%rename  (__iless__) operator-=;
%rename  (__imul__) operator*=;
%rename  (__idiv__) operator/=;
%rename  (__eq__) operator==;
%rename  (__ne__) operator!=;


%rename("$ignore", regexmatch$name="^createSPtr.*") "";

// swig see those classes as abstract => no wrappers for constructors
// note: this means something may be wrong in .hpp headers => use -Wall to
// detect it

// swig does not see the using Simulation::update ?
%feature("notabstract") TimeStepping;
%feature("notabstract") EventDriven;




// includes and imports
%include "SiconosPointers.hpp"
%include "SiconosVisitables.hpp"

%import "SiconosVisitor.hpp"
%import "SiconosSerialization.hpp"

%import "SiconosProperties.hpp"

 /* swig has difficulties with this macro in SiconosProperties */
#undef INSTALL_GRAPH_PROPERTIES
#define INSTALL_GRAPH_PROPERTIES(X,Y)


%include "SiconosGraph.hpp"


TYPEDEF_SPTR(_DynamicalSystemsGraph);
%feature("director") _DynamicalSystemsGraph;
%shared_ptr( SiconosGraph<std11::shared_ptr<DynamicalSystem>, 
                          std11::shared_ptr<Interaction>, 
                          SystemProperties, InteractionProperties, 
                          GraphProperties >);

TYPEDEF_SPTR(_InteractionsGraph);
%feature("director") _InteractionsGraph;
%shared_ptr( SiconosGraph<std11::shared_ptr<Interaction>, 
                          std11::shared_ptr<DynamicalSystem>, 
                          InteractionProperties, SystemProperties, 
                          GraphProperties >);

TYPEDEF_SPTR(DynamicalSystemsGraph);
%feature("director") DynamicalSystemsGraph;
%shared_ptr(DynamicalSystemsGraph);

TYPEDEF_SPTR(InteractionsGraph);
%feature("director") InteractionsGraph;
%shared_ptr(InteractionsGraph);

// must be specified after %shared_ptr, if ever needed
%template(_DynamicalSystemsGraph) SiconosGraph<
  std11::shared_ptr<DynamicalSystem>, 
  std11::shared_ptr<Interaction>, 
  SystemProperties, InteractionProperties, 
  GraphProperties >;

%template(SP_DynamicalSystemsGraph) std11::shared_ptr<
  SiconosGraph<std11::shared_ptr<DynamicalSystem>, 
               std11::shared_ptr<Interaction>, 
               SystemProperties, InteractionProperties, 
               GraphProperties > >;

%template(_InteractionsGraph) SiconosGraph<
  std11::shared_ptr<Interaction>, 
  std11::shared_ptr<DynamicalSystem>, 
  InteractionProperties, SystemProperties, 
  GraphProperties >;

%template(SP_InteractionsGraph) std11::shared_ptr<
  SiconosGraph<std11::shared_ptr<Interaction>, 
               std11::shared_ptr<DynamicalSystem>, 
               InteractionProperties, SystemProperties, 
               GraphProperties > >;

%rename (ioMatrix_read) ioMatrix::read; 
%rename (ioMatrix_write) ioMatrix::write; 
%include "ioMatrix.hpp"

%include "SimulationTypeDef.hpp" 

%include "InteractionsSet.hpp"
%include "SiconosSet.hpp"

%import "boost/config.hpp"
%import "boost/graph/graph_utility.hpp"

%include "Tools.hpp"

%include "addons.hpp"

%include "KernelRegistration.i"


%define PY_REGISTER(TYPE)
%rename  (__getitem__) TYPE ## ::operator[];
%rename  (__add__) TYPE ## ::operator+;
%rename  (__mul__) TYPE ## ::operator*;
%rename  (__div__) TYPE ## ::operator/;
%rename  (__iadd__) TYPE ## ::operator+=;
%rename  (__imul__) TYPE ## ::operator*=;
%rename  (__idiv__) TYPE ## ::operator/=;
%rename  (__eq__) TYPE ## ::operator==;
%rename  (__ne__) TYPE ## ::operator!=;
%rename  (__copy__) TYPE ## ::operator=;
%feature("director") TYPE;
%ignore STD11::enable_shared_from_this<TYPE>;
%template (shared ## TYPE) STD11::enable_shared_from_this<TYPE>;
%typemap(directorin) TYPE& ()
{
  // %typemap(directorin) (TYPE&) ()
  // swig issue shared pointer check in wrappers even if arg is a ref
  SP::TYPE myptemp(createSPtr##TYPE($1));
  $input = SWIG_NewPointerObj(SWIG_as_voidptr(&myptemp), 
                              SWIGTYPE_p_std11__shared_ptrT_##TYPE##_t, 0);
}
%shared_ptr(TYPE); 
%enddef


%define PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%rename  (__getitem__) TYPE ## ::operator[];
%rename  (__add__) TYPE ## ::operator+;
%rename  (__mul__) TYPE ## ::operator*;
%rename  (__div__) TYPE ## ::operator/;
%rename  (__iadd__) TYPE ## ::operator+=;
%rename  (__imul__) TYPE ## ::operator*=;
%rename  (__idiv__) TYPE ## ::operator/=;
%rename  (__eq__) TYPE ## ::operator==;
%rename  (__ne__) TYPE ## ::operator!=;
%rename  (__copy__) TYPE ## ::operator=;
%ignore STD11::enable_shared_from_this<TYPE>;
%template (shared ## TYPE) STD11::enable_shared_from_this<TYPE>;
%shared_ptr(TYPE); 
%enddef

%define PY_REGISTER_BULLET_COLLISION_DETECTION(X) 
TYPEDEF_SPTR(X);
PY_REGISTER_WITHOUT_DIRECTOR(X)
%enddef

%define PY_REGISTER_BULLET_LINEAR_MATH(X)
TYPEDEF_SPTR(X);
PY_REGISTER_WITHOUT_DIRECTOR(X)
%enddef

 // registered classes in KernelRegistration.i
KERNEL_REGISTRATION();

// ignores

// Bullet
// (mostly because not defined in <name>.h)
%ignore btCapsuleShapeX;
%ignore btCapsuleShapeZ;
%ignore btConeShapeX;
%ignore btConeShapeZ;
%ignore btCylinderShapeX;
%ignore btCylinderShapeZ;
%ignore btConvexInternalAabbCachingShape;
%ignore btPolyhedralConvexAabbCachingShape;
%ignore btBU_Simplex1to4;
%ignore m_vertices1;

%ignore btVector4;



// createSPtr*
%ignore nullDeleter;

// defined in SiconosVector.cpp
%ignore setBlock;
%ignore add;
%ignore sub;
%ignore axpby;
%ignore axpy;
%ignore inner_prod;
%ignore outer_prod;
%ignore scal;
%ignore subscal;
%ignore cross_product;

// defined in SimpleMatrix.cpp 
%ignore private_addprod;
%ignore private_prod;
%ignore prod;
%ignore axpy_prod;
%ignore subprod;
%ignore axpy_prod;
%ignore gemv;
%ignore gemm;

%ignore getWMap;
%ignore getWBoundaryConditionsMap;
%ignore getDSBlocks;
%ignore getInvMSimple;
%ignore getInvMBlock;

// do not wrap visitor visit : this lead to a huge amount of wrapper
// code generation and this fail at compile time on shared_ptr freearg
%ignore SiconosVisitor::visit;

// cannot compile wrapper
%ignore statOut;

%include "SiconosAlgebraTypeDef.hpp"
%include "SiconosAlgebra.hpp"

%import "RelationNamespace.hpp";

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



%inline
%{
  
  /* Note: without the PyCObject stuff the python
   * wrapper fail on this, the numpy vector points on a deleted
   * memory!*/
  
  const SP::SiconosVector getVector(SP::SiconosVector v)
  {
    return v;
  };

  const SP::SiconosMatrix getMatrix(SP::SiconosMatrix v)
  {
    return v;
  };
  
  /* to make swig define SWIGTYPE_p_PyArrayObject */
  const PyArrayObject* getVector(PyArrayObject* v)
  {
    return v;
  };


%}




// needed templates

%template (InteractionsSet) SiconosSet<Interaction,double*>;

%template (dsi) std::pair<unsigned int, unsigned int >;

%template (dsp) std::pair<std11::shared_ptr<DynamicalSystem>, 
                          std11::shared_ptr<DynamicalSystem> >;

%template (dspv) std::vector<std::pair<std11::shared_ptr<DynamicalSystem>, 
                                       std11::shared_ptr<DynamicalSystem> > >;

%template (dsiv) std::vector<std::pair<unsigned int, unsigned int > >;

%template(unsignedintv) std11::shared_ptr<std::vector<unsigned int> >;

// not sufficient
%ignore Question<bool>;
%template (qbool) Question<bool>;

%ignore Question<unsigned int>;
%template (quint) Question<unsigned int>;

// suppress warning
%ignore  STD11::enable_shared_from_this< Hashed >;
%template (sharedHashed) STD11::enable_shared_from_this< Hashed >;


// include registered headers

#undef PY_REGISTER
%define PY_REGISTER(X)
%include "X.hpp";
%enddef


#undef PY_REGISTER_WITHOUT_DIRECTOR
%define PY_REGISTER_WITHOUT_DIRECTOR(X)
%include "X.hpp";
%enddef

#undef PY_REGISTER_BULLET_COLLISION_DETECTION
%define PY_REGISTER_BULLET_COLLISION_DETECTION(X)
%include "BulletCollision/CollisionShapes/X.h";
%enddef

#undef PY_REGISTER_BULLET_LINEAR_MATH
%define PY_REGISTER_BULLET_LINEAR_MATH(X)
%include "LinearMath/X.h";
%enddef

%shared_ptr(_SolverOptions);
TYPEDEF_SPTR(_SolverOptions);

// note : deleteSolverOptions is call by ~LCP(), ~FrictionContact(), etc.

%include "SolverOptions.h"

KERNEL_REGISTRATION();

%include "FrontEndConfig.h";

#ifdef HAVE_SICONOS_IO
%include "SiconosRestart.hpp";
#endif

%fragment("StdSequenceTraits");

%fragment("StdMapTraits");

