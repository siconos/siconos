// -*- C++ -*-
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.
//
// Copyright 2016 INRIA.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

// SWIG interface for Siconos Kernel
%module(package="siconos", directors="1", allprotected="1") kernel

%include start.i

// generated docstrings from doxygen xml output
%include kernel-docstrings.i

#ifdef WITH_SERIALIZATION
%{
#define KERNEL_ONLY
#include <SiconosFull.hpp>
%}
#endif
%include picklable.i

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
// remove the visitor stuff, I see no use of it in Python for now --xhub
%ignore *::acceptSerializer;
%ignore *::acceptType;
%ignore *::accept;
%ignore *::acceptSP;
// do not wrap visitor visit : this lead to a huge amount of wrapper
// code generation and this fail at compile time on shared_ptr freearg
%ignore SiconosVisitor;

%ignore visit;

%ignore Type::str;

// cannot compile wrapper
%ignore statOut;

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
%ignore getInvMBlock;

%warnfilter(509) rotateAbsToBody;
%warnfilter(509) changeFrameAbsToBody;
%warnfilter(509) changeFrameBodyToAbs;


 // common declarations with upper modules : Mechanics, IO, ...
%include handleException.i

%include sharedPointers.i

// handle stl data types
%include stl.i

// 1. Vector and Matrix <=> numpy array (dense only)
%include SiconosAlgebra.i


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



// swig see those classes as abstract => no wrappers for constructors
// note: this means something may be wrong in .hpp headers => use -Wall to
// detect it

// yes, undefined private copy constructors
%feature("notabstract") TimeStepping;
%feature("notabstract") TimeSteppingCombinedProjection;
%feature("notabstract") TimeSteppingDirectProjection;
%feature("notabstract") EventDriven;

// common declarations with Numerics

// note : solver_options_delete is call by ~LCP(), ~FrictionContact(), etc.
%shared_ptr(SolverOptions);
%shared_ptr(NumericsMatrix);
%shared_ptr(CSparseMatrix);
%shared_ptr(SparseBlockStructuredMatrix);
%shared_ptr(GlobalFrictionContactProblem);

%include solverOptions.i

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
      $result = SWIG_NewPointerObj(SWIG_as_voidptr($1->matrix1), $descriptor(SparseBlockStructuredMatrix *), 0);
    }
    else // failing silently does not seem a good option
    {
      Py_INCREF(Py_None);
      $result = Py_None;
    }
  }
}
%import NumericsMatrix.h
%import SparseMatrix.h
%import SparseBlockMatrix.h

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

%import "boost/config.hpp"
%import "boost/graph/graph_utility.hpp"

%include "Tools.hpp"

%include "addons.hpp"

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

  SP::RelayNSL cast_RelayNSL(SP::NonSmoothLaw nslaw)
  {
    return std11::dynamic_pointer_cast<RelayNSL>(nslaw);
  }

  SP::NewtonImpactNSL cast_NewtonImpactNSL(SP::NonSmoothLaw nslaw)
  {
    return std11::dynamic_pointer_cast<NewtonImpactNSL>(nslaw);
  }

%}

//namespace std {

  %template (dspv) std::vector<std::pair<std11::shared_ptr<DynamicalSystem>,
                                         std11::shared_ptr<DynamicalSystem> > >;

  %template (dsiv) std::vector<std::pair<unsigned int, unsigned int > >;


  %template (dsi) std::pair<unsigned int, unsigned int >;

  %template (dsp) std::pair<std11::shared_ptr<DynamicalSystem>,
                            std11::shared_ptr<DynamicalSystem> >;

//BouncingBallNETS.py, attempt to reach DSlink as a vector...
//swig failure.
//%shared_ptr(VectorOfBlockVectors);
//%template (vectorOfBlockVectors) std::vector<std11::shared_ptr<BlockVector> >;
///
//}


%template(unsignedintv) std11::shared_ptr<std::vector<unsigned int> >;

// not sufficient
%ignore Question<bool>;
%template (qbool) Question<bool>;

%ignore Question<unsigned int>;
%template (quint) Question<unsigned int>;


%ignore OSNSMatrix::updateSizeAndPositions;

// registered classes in KernelRegistration.i

%include KernelRegistration.i
%include pyRegister.i
KERNEL_REGISTRATION();

%include pyInclude.i

KERNEL_REGISTRATION()

%fragment("StdSequenceTraits");

%fragment("StdMapTraits");


