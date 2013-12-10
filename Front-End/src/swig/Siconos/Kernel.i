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

 // signatures
%feature("autodoc", 1);

// generated docstrings from doxygen xml output
%include Kernel-docstrings.i

%include start.i

%{
#include <SiconosKernel.hpp>
#include <SiconosAlgebra.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <SiconosVisitor.hpp>
#include "addons.hpp"
#include <boost/type_traits/is_polymorphic.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/typeof/typeof.hpp>
%}

// ignores
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

%ignore visit;

%ignore Type::str;

// cannot compile wrapper
%ignore statOut;

// mandatory !
%rename (lambda_) lambda;

 // common declarations with upper modules : Mechanics, IO, ...
%include handleException.i

%include sharedPointers.i

// handle stl data types
%include stl.i

// 1. Vector and Matrix <=> numpy array (dense only)
%include KernelTypes.i


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


// swig see those classes as abstract => no wrappers for constructors
// note: this means something may be wrong in .hpp headers => use -Wall to
// detect it

// yes, undefined private copy constructors
%feature("notabstract") TimeStepping;
%feature("notabstract") EventDriven;

// common declarations with Numerics

// note : deleteSolverOptions is call by ~LCP(), ~FrictionContact(), etc.
%shared_ptr(_SolverOptions);
%shared_ptr(NumericsOptions);
%shared_ptr(NumericsMatrix);
%shared_ptr(SparseMatrix);
%shared_ptr(SparseBlockStructuredMatrix);
%shared_ptr(FrictionContactProblem);
%shared_ptr(GlobalFrictionContactProblem);

%include NumericsOptions.h
%include solverOptions.i
%include FrictionContactProblem.i

// access NumericsMatrix cf Numerics.i
%typemap(out) (std11::shared_ptr<NumericsMatrix>) {
  npy_intp dims[2];

  if (!$1) 
  {
    Py_INCREF(Py_None);
    $result = Py_None;
  }
  else
  {
    dims[0] = $1->size0;
    dims[1] = $1->size1;
    
    if ($1->matrix0)
    {
      PyObject *obj = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, $1->matrix0);
      PyArrayObject *array = (PyArrayObject*) obj;
      if (!array || !require_fortran(array)) SWIG_fail;
      $result = obj;
    }
    else if($1->matrix1)
    {
      // matrix is sparse : return opaque pointer
      $result = SWIG_NewPointerObj(SWIG_as_voidptr($1->matrix1), SWIGTYPE_p_SparseBlockStructuredMatrix, 0);
    }
    else // failing silently does not seem a good option
    {
      Py_INCREF(Py_None);
      $result = Py_None;
    }
  }
}
%include NumericsMatrix.h
%include SparseMatrix.h
%include SparseBlockMatrix.h

 // segfaults...
 // we cannot share data struct
 //%import Numerics.i

%include "SiconosConst.hpp"

%include "SiconosVisitables.hpp"

%import "SiconosVisitor.hpp"
%import "Question.hpp"
%import "TypeName.hpp"

%import "SiconosSerialization.hpp"

%import "SiconosProperties.hpp"

%include graph.i

%rename (ioMatrix_read) ioMatrix::read; 
%rename (ioMatrix_write) ioMatrix::write; 
%include "ioMatrix.hpp"

%include "SimulationTypeDef.hpp" 

%include "SiconosSet.hpp"

%import "boost/config.hpp"
%import "boost/graph/graph_utility.hpp"

%include "ControlTypeDef.hpp"

%include "Tools.hpp"

%include "addons.hpp"

%include "ControlTools.hpp"

// fix : how to prevent swig to generate getter/setter for mpz_t ?
// on some distrib %import gmp.h is not sufficient as gmp-<arch>.h may be used
typedef struct
{} __mpz_struct;
typedef __mpz_struct mpz_t[1];

%include "SiconosAlgebraTypeDef.hpp"
%include "SiconosAlgebra.hpp"

%import "RelationNamespace.hpp";




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

  SP::NewtonImpactFrictionNSL cast_NewtonImpactFrictionNSL(SP::NonSmoothLaw nslaw)
  {
    return std11::dynamic_pointer_cast<NewtonImpactFrictionNSL>(nslaw);
  }

%}

//namespace std {

  %template (dspv) std::vector<std::pair<std11::shared_ptr<DynamicalSystem>, 
                                         std11::shared_ptr<DynamicalSystem> > >;

  %template (dsiv) std::vector<std::pair<unsigned int, unsigned int > >;


  %template (dsi) std::pair<unsigned int, unsigned int >;

  %template (dsp) std::pair<std11::shared_ptr<DynamicalSystem>, 
                            std11::shared_ptr<DynamicalSystem> >;

//}


%template(unsignedintv) std11::shared_ptr<std::vector<unsigned int> >;

// not sufficient
%ignore Question<bool>;
%template (qbool) Question<bool>;

%ignore Question<unsigned int>;
%template (quint) Question<unsigned int>;

// suppress warning
%ignore  STD11::enable_shared_from_this< Hashed >;
%template (sharedHashed) STD11::enable_shared_from_this< Hashed >;



%ignore OSNSMatrix::updateSizeAndPositions;

// registered classes in KernelRegistration.i

%include KernelRegistration.i
%include pyRegister.i
KERNEL_REGISTRATION();

%include pyInclude.i

KERNEL_REGISTRATION()

%fragment("StdSequenceTraits");

%fragment("StdMapTraits");

