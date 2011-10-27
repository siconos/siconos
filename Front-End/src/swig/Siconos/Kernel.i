// -*- C++ -*-
// Siconos-Front-End , Copyright INRIA 2005-2011.
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
%module(directors="1") Kernel

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
#include "SiconosKernel.hpp"
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

#define DEBUG_MESSAGES 1
#include <debug.h>

%} 

// mandatory !
%rename (lambda_) lambda;

// shared ptr management
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
  static int sequenceToUnsignedIntVector(PyObject *input, boost::shared_ptr<std::vector<unsigned int> > ptr) 
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

// boost >= 1.40
%import "boost/version.hpp"
#if (BOOST_VERSION >= 104000)
%ignore boost::enable_shared_from_this::operator=;
%import "boost/smart_ptr/enable_shared_from_this.hpp"
#else
%import "boost/enable_shared_from_this.hpp"
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

%include "SiconosPointers.hpp"
%include "SiconosVisitables.hpp"
%import "SiconosVisitor.hpp"
%import "SiconosSerialization.hpp"
%import "ioObject.hpp"
%include "ioMatrix.hpp"

%import "SiconosGraph.hpp"

%include "SimulationTypeDef.hpp" 

%include "InteractionsSet.hpp"
%include "SiconosSet.hpp"

%import "boost/config.hpp"
%import "boost/graph/graph_utility.hpp"

// what's wrong with this ignore ??
%ignore PURE_DEF;
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
%ignore boost::enable_shared_from_this<TYPE>;
%template (shared ## TYPE) boost::enable_shared_from_this<TYPE>;
%shared_ptr(TYPE); 
%enddef


 // registered classes in KernelRegistration.i
KERNEL_REGISTRATION();

// ignores

%ignore nullDeleter;

// defined in SimpleVector.cpp
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
    boost::shared_ptr<void> ref;

    SharedPointerKeeper(boost::shared_ptr<void> v) : ref(v) 
    {
      DEBUG_PRINTF("SharedPointerKeeper : use_count %ld",v.use_count());
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
       new SharedPointerKeeper(boost::static_pointer_cast<void>(mysharedptr));
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

  
  /* for testing purpose : without the PyCObject stuff the python
   * wrapper fail on this, the numpy vector points on a deleted
   * memory!*/
  const SP::SimpleVector getVector(SP::SimpleVector v)
  {
    return v;
  };
  
  /* to make swig define SWIGTYPE_p_PyArrayObject */
  const PyArrayObject* getVector(PyArrayObject* v)
  {
    return v;
  };
  

%}


// include registered headers

#undef PY_REGISTER
%define PY_REGISTER(X)
%include "X.hpp";
%enddef


%template (InteractionsSet) SiconosSet<Interaction,double*>;

%template (dsi) std::pair<unsigned int, unsigned int >;

%template (dsp) std::pair<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<DynamicalSystem> >;

%template (dspv) std::vector<std::pair<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<DynamicalSystem> > >;

%template (dsiv) std::vector<std::pair<unsigned int, unsigned int > >;

%template (ioMatrix) ioObject<SiconosMatrix>; 


// not sufficient
%ignore Question<bool>;
%template (qbool) Question<bool>;

%ignore Question<unsigned int>;
%template (quint) Question<unsigned int>;

// suppress warning
%ignore  boost::enable_shared_from_this< Hashed >;
%template (sharedHashed) boost::enable_shared_from_this< Hashed >;


KERNEL_REGISTRATION();


%fragment("StdSequenceTraits");

%fragment("StdMapTraits");

