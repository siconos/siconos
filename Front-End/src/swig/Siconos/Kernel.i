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

// Siconos.i - SWIG interface for Siconos
%module Kernel

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
%} 

// shared ptr management
%include "boost_shared_ptr.i"

// numpy macros
%include numpy.i 	

%init %{
  import_array();
%}

// Handle standard exceptions
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
}

// 1. Vector and Matrix <=> numpy array (dense only)


// numpy array to SP::SimpleVector (here a SiconosVector is always a
// SimpleVector)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosVector>)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in,fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> (PyArrayObject* array=NULL, int is_new_object)
{
  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,1) ||
      !require_native(array) || !require_contiguous(array)) SWIG_fail;

  SP::SimpleVector tmp;
  tmp.reset(new SimpleVector(array_size(array,0)));
  // copy : with SimpleVector based on resizable std::vector there is
  // no other way
  memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
  $1 = tmp;
 }

// numpy array to SP::SimpleMatrix (here a SiconosMatrix is always a
// SimpleMatrix)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosMatrix>)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in) boost::shared_ptr<SiconosMatrix> (PyArrayObject* array=NULL, int is_new_object) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,2) ||
      !require_native(array) || !require_contiguous(array)) SWIG_fail;

  SP::SimpleMatrix tmp;
  tmp.reset(new SimpleMatrix(array_size(array,0), array_size(array,1)));
  // copy : with SimpleMatrix based on resizable std::vector
  memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
  $1 = tmp;
 }


// from C++ to python 
%template() boost::shared_ptr<SiconosVector>;
%template() boost::shared_ptr<SimpleVector>;
%template() boost::shared_ptr<SiconosMatrix>;
%template() boost::shared_ptr<SimpleMatrix>;

%typemap(out) boost::shared_ptr<SiconosVector>
{
  int this_vector_dim[1];
  this_vector_dim[0]=$1->size();
  $result = PyArray_SimpleNewFromData(1,this_vector_dim,NPY_DOUBLE,$1->getArray());
}

%typemap(out) boost::shared_ptr<SiconosMatrix>
{
  int this_matrix_dim[2];
  this_matrix_dim[0]=$1->size(0);
  this_matrix_dim[1]=$1->size(1);
  $result = PyArray_SimpleNewFromData(2,this_matrix_dim,NPY_DOUBLE,$1->getArray());
  PyArray_UpdateFlags((PyArrayObject *)$result, NPY_FORTRAN);
}



// SiconosMatrix in api mean SimpleMatrix here
%apply (boost::shared_ptr<SiconosVector>) { (SP::SimpleVector) };
%apply (boost::shared_ptr<SiconosVector>) { (boost::shared_ptr<SimpleVector>) };
%apply (boost::shared_ptr<SiconosVector>) { (SP::SiconosVector) };

%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SimpleMatrix) };
%apply (boost::shared_ptr<SiconosMatrix>) { (boost::shared_ptr<SimpleMatrix>) };
%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SiconosMatrix) }; 


// 2. try to hide SP::Type on python side

%define SP_TYPE(TYPE)
%typecheck(SWIG_TYPECHECK_POINTER)
(SP::##TYPE)
{
  void *p;
  int r=0;
  r=SWIG_ConvertPtr($input, &p, SWIGTYPE_p_##TYPE, 0 | 0);
  if(!SWIG_IsOK(r))
  {
    r=SWIG_ConvertPtr($input, &p, SWIGTYPE_p_boost__shared_ptrT_##TYPE##_t, 0 | 0);
    if(!SWIG_IsOK(r)) SWIG_exception_fail(SWIG_ArgError(r), "not a TYPE or SP:: TYPE");
  }
};

%typemap(in) SP::##TYPE (PyObject* obj, void *p, int r)
{
  r = SWIG_ConvertPtr($input, &p, SWIGTYPE_p_##TYPE, 0 | 0);
  if(!SWIG_IsOK(r))
  {
    r=SWIG_ConvertPtr($input, &p, SWIGTYPE_p_boost__shared_ptrT_##TYPE##_t, 0 | 0);
    if(!SWIG_IsOK(r)) SWIG_exception_fail(SWIG_ArgError(r), "not a TYPE or SP:: TYPE");
    SP::##TYPE * ptmp = reinterpret_cast< SP::##TYPE * >(p);
    $1 = *ptmp;
  }
  else
  {
    TYPE *pr;
    pr = reinterpret_cast< TYPE * > (p);
    
    SP::##TYPE tmp = createSPtr##TYPE(*pr);
    $1 = tmp;
    
  }
};
%enddef


SP_TYPE(NonSmoothLaw);
SP_TYPE(Relation);


// dummy namespaces to make swig happy
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


namespace SP
{
};


// note Visitor etc. => module Visitor, Type if ever needed

%ignore PURE_DEF;
%ignore VIRTUAL_ACCEPT_VISITORS;
%ignore ACCEPT_STD_VISITORS;
%ignore ACCEPT_VISITORS;

%ignore nullDeleter;

%ignore DynamicalSystem;

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
%ignore operator +;
%ignore operator [];
%ignore operator -;
%ignore operator *;
%ignore operator /;
%ignore operator ==;
%ignore operator =;

// defined in SimpleMatrix.cpp 
%ignore private_addprod;
%ignore private_prod;
%ignore prod;
%ignore axpy_prod;
%ignore subprod;
%ignore axpy_prod;
%ignore gemv;
%ignore gemm;



%include "SiconosPointers.hpp"
%include "SiconosGraph.hpp"


DEFINE_SPTR(SiconosVector);
DEFINE_SPTR(SiconosMatrix);

DEFINE_SPTR(NonSmoothLaw);
DEFINE_SPTR(Interaction);
DEFINE_SPTR(Relation);

%include "SimpleVector.hpp"
%include "SimpleMatrix.hpp"
%include "DynamicalSystem.hpp"
%include "LagrangianDS.hpp"
%include "LagrangianLinearTIDS.hpp"
%include "NewtonImpactNSL.hpp"
%include "Relation.hpp"
%include "LagrangianR.hpp"
%include "LagrangianLinearTIR.hpp"
%include "Interaction.hpp"
%include "Model.hpp"
%include "OneStepIntegrator.hpp"
%include "Moreau.hpp"
%include "LCP.hpp"
%include "TimeStepping.hpp"

