// -*- C++ -*-
// Siconos-Front-End version 3.2.0, Copyright INRIA 2005-2010.
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


 // where is doxygen feature ?
%feature("autodoc", "1");

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
%include "exception.i"
%exception
{
  try
  {
    $action
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch (const std::out_of_range& e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
  catch (const SiconosException& e)
  {
    SWIG_exception(SWIG_SystemError, e.report().c_str());
  }
  catch (const Swig::DirectorException& e)
  {
    SWIG_exception(SWIG_ValueError, e.getMessage());
  }
}

// handle stl data types
%include "stl.i"

// 1. Vector and Matrix <=> numpy array (dense only)


%include "KernelTypes.i"

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


%include "SiconosVisitables.hpp"
%import "SiconosVisitor.hpp"

%include "SiconosPointers.hpp"
%import "SimulationTypeDef.hpp" 

%include "InteractionsSet.hpp"
%include "SiconosSet.hpp"

%import "SiconosGraph.hpp"


%import "boost/config.hpp"
%import "boost/graph/graph_utility.hpp"


%include "Tools.hpp"
%include "addons.hpp"

%include "KernelRegistration.i"

#define PY_REGISTER(TYPE) \
%rename  (__getitem__) TYPE ## ::operator[]; \
%rename  (__add__) TYPE ## ::operator+; \
%rename  (__mul__) TYPE ## ::operator*; \
%rename  (__div__) TYPE ## ::operator/; \
%rename  (__iadd__) TYPE ## ::operator+=; \
%rename  (__imul__) TYPE ## ::operator*=; \
%rename  (__idiv__) TYPE ## ::operator/=; \
%rename  (__eq__) TYPE ## ::operator==; \
%rename  (__ne__) TYPE ## ::operator!=; \
%rename  (__copy__) TYPE ## ::operator=; \
%feature("director") TYPE; \
%shared_ptr(TYPE); 


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



// include registered headers

#undef PY_REGISTER
%define PY_REGISTER(X)
%include "X.hpp";
%enddef

KERNEL_REGISTRATION();





%fragment("StdSequenceTraits");

%fragment("StdMapTraits");


%template (InteractionsSet) SiconosSet<Interaction,double*>;

%template (sharedModel) boost::enable_shared_from_this<Model>;

%template (dsv) std::vector<boost::shared_ptr<DynamicalSystem> >;

%template (dsi) std::pair<unsigned int, unsigned int >;

%template (dsp) std::pair<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<DynamicalSystem> >;

%template (dspv) std::vector<std::pair<boost::shared_ptr<DynamicalSystem>, boost::shared_ptr<DynamicalSystem> > >;

%template (dsiv) std::vector<std::pair<unsigned int, unsigned int > >;


