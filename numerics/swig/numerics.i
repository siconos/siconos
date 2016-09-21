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

// numerics.i - SWIG interface for siconos numerics component.
%module(package="siconos") numerics

// basics, mostly numpy.i stuff. 
%include start.i

// generated docstrings from doxygen xml output
%include numerics-docstrings.i

// list of functions that returns objects that they have created
%include numerics_newobjects.i

%include <typemaps/swigmacros.swg>

%define %SN_INPUT_CHECK_RETURN(var, output, type_)
{ void* _ptr = NULL; if (!SWIG_IsOK(SWIG_ConvertPtr(var, &_ptr, SWIGTYPE_p_##type_, 0 |  0 ))) {
    char errmsg[1024];
    snprintf(errmsg, sizeof(errmsg), "Argument check failed! Argument %s has the wrong type, should be %s", #var, #type_);
    PyErr_SetString(PyExc_TypeError, errmsg); PyErr_PrintEx(0); return NULL; }
    output = %static_cast(_ptr, type_*);
   }
%enddef

%{
#include "SiconosNumerics.h"
#include "SiconosConfig.h"
#include "SolverOptions.h"
#include "SparseMatrix.h"
#include "NumericsMatrix.h"
#include "SparseBlockMatrix.h"
#include "NumericsSparseMatrix.h"
#include "Numerics_functions.h"
#include "SiconosSets.h"
#include "GAMSlink.h"
#include "NumericsFwd.h"
  %}

#ifdef WITH_SERIALIZATION
%{
#include <SiconosFullNumerics.hpp>
%}
#endif
%include picklable.i

// needed macros
%include "SiconosConfig.h"

// declare C++ shared_ptrs to C structs
// need to do this here for other modules to reference numerics structs by shared_ptr.
// swig requires same namespace 'std11' is used.
%{
#if defined(SICONOS_STD_SHARED_PTR) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
namespace std11 = std;
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
namespace std11 = boost;
#endif
%}
#define SWIG_SHARED_PTR_NAMESPACE std11
%include boost_shared_ptr.i

// put needed shared_ptrs here
// commented-out for now, swig insists on calling freeFrictionContactProblem()
// instead of respecting the shared_ptr!
// %shared_ptr(FrictionContactProblem)


 // more convenient
 %rename (LCP) LinearComplementarityProblem;
 %rename (MLCP) MixedLinearComplementarityProblem;
 %rename (MCP) MixedComplementarityProblem;
 %rename (VI) VariationalInequality;
 %rename (AVI) AffineVariationalInequalities;

 %ignore lcp_compute_error_only;

 // -- Numpy typemaps --
 // See http://docs.scipy.org/doc/numpy/reference/swig.interface-file.html.
 %apply (int DIM1 , int DIM2 , double* IN_FARRAY2)
 {(int nrows, int ncols, double* data          )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizex, double *x          )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizey, double *y          )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizeq, double *q         )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizez, double *z        )};

 %apply (int DIM1 , double* INPLACE_ARRAY1)
 {(int sizew, double *w        )};

 %apply (double IN_ARRAY1[ANY])
 {(double reaction[3])}

 %apply (double IN_ARRAY1[ANY])
 {(double velocity[3])}

 %apply (double IN_ARRAY1[ANY])
 {(double rho[3])}

 %apply (double ARGOUT_ARRAY1[ANY])
 {(double ACresult[3])}

 %apply (double ARGOUT_ARRAY1[ANY])
 {(double result[3])}

 %apply (double ARGOUT_ARRAY1[ANY])
 {(double A[9])}

 %apply (double ARGOUT_ARRAY1[ANY])
 {(double B[9])}


%include Numerics_typemaps_problems.i
%include Numerics_typemaps_basic_linalg.i
%include Numerics_typemaps_numericsmatrices.i
%include NonSmoothDrivers.h
%include solverOptions.i

// this has to die --xhub
// info param
%typemap(in, numinputs=0) (int *info) (int temp_info = -1)
{
  // a default initialization : solver may stop if *info = 0 (checkTrivialCase)
  // checkTrivialCase => better if directly in solvers, not in driver.
  $1 = &temp_info;
}

%inline
%{
 static unsigned int isqrt(unsigned int n)
  {
    unsigned int c = 0x8000;
    unsigned int g = 0x8000;

    for(;;)
    {
      if(g*g > n)
      g ^= c;
      c >>= 1;
      if(c == 0)
        return g;
      g |= c;
    }
  }

  static int compiled_in_debug_mode()
  {
#ifdef NDEBUG
    return 0;
#else
    return 1;
#endif
  }

#ifdef __cplusplus
  extern "C" {
#endif
  // this can't be static --xhub
  void set_cstruct(uintptr_t p_env, void* p_struct);
  void set_cstruct(uintptr_t p_env, void* p_struct)
  {
    *(void**)p_env = p_struct;
  }
#ifdef __cplusplus
}
#endif
%}

%fragment("NumPy_Fragments");

// includes in 'begin' mandatory to avoid mess with
// solverOptions.i, numerics_common and fwd decl
// all this because of SolverOptions extend.
%begin %{
#include "relay_cst.h"
#include "AVI_cst.h"
#include "SOCLCP_cst.h"
#include "Friction_cst.h"
#include "lcp_cst.h"
#include "MCP_cst.h"
#include "mlcp_cst.h"
#include "NCP_cst.h"
#include "VI_cst.h"
#include "GenericMechanical_cst.h"
#include "fc2d_Solvers.h"
#include "fc3d_Solvers.h"
#include "gfc3d_Solvers.h"
#include "MCP_Solvers.h"
#include "MLCP_Solvers.h"
#include "NonSmoothDrivers.h"  
  %}

//Relay
%include "relay_cst.h"

%include numerics_MLCP.i

/////////////////////////
// This is common to all the problem defined below
// MLCP should be defined above this statement
/////////////////////////

%typemap(out) (double* q) {
  npy_intp dims[1];

  if (!arg1->M) { PyErr_SetString(PyExc_TypeError, "M is not present, don't known the size"); SWIG_fail; }

  dims[0] = arg1->M->size0;
  if ($1)
  {
    PyObject *obj = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, $1);
    PyArrayObject *array = (PyArrayObject*) obj;
    if (!array || !require_fortran(array)) SWIG_fail;
    $result = obj;
  }
  else
    SWIG_fail;
 }

%typemap(out) (double* mu) {
  npy_intp dims[1];

  if (arg1->numberOfContacts <= 0) { PyErr_SetString(PyExc_TypeError, "numberOfContacts is not set"); SWIG_fail; }

  dims[0] = arg1->numberOfContacts;
  if ($1)
  {
    PyObject *obj = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, $1);
    PyArrayObject *array = (PyArrayObject*) obj;
    if (!array || !require_fortran(array)) SWIG_fail;
    $result = obj;
  }
  else
    SWIG_fail;
 }

// vectors of size numberOfContacts
%typemap(memberin) (double *mu) {
  // Still some dark magic :( --xhub
  if (arg1->numberOfContacts <= 0)
  {
    PyErr_SetString(PyExc_RuntimeError, "numberOfContacts is not set, it sould be done first!");
    SWIG_fail;
  }

  if (arg1->numberOfContacts !=  array_size(array2, 0))
  {
    char msg[1024];
    snprintf(msg, sizeof(msg), "Size of mu is %ld, but the number of contacts is %d! Both should be equal!\n", array_size(array2, 0), arg1->numberOfContacts);
    PyErr_SetString(PyExc_RuntimeError, msg);
    SWIG_fail;
  }

  if (!$1) { $1 = (double*)malloc(arg1->numberOfContacts * sizeof(double)); }
  memcpy($1, $input, arg1->numberOfContacts * sizeof(double));

 }

// vectors of size M
%typemap(memberin) (double *q) {
  // Still some dark magic :( --xhub
 char msg[1024];
  assert(arg1);
  if (!arg1->M)
  {
    PyErr_SetString(PyExc_RuntimeError, "M is not initialized, it sould be done first!");
    SWIG_fail;
  }

  int size = arg1->M->size0;
  if (size !=  array_size(array2, 0))
  {
    snprintf(msg, sizeof(msg), "Size of q is %ld, but the size of M is %d! Both should be equal!\n", array_size(array2, 0), size);
    PyErr_SetString(PyExc_RuntimeError, msg);
    SWIG_fail;
  }

  if (!$1) { $1 = (double*)malloc(size * sizeof(double)); }
  memcpy($1, $input, size * sizeof(double));

 }


%include numerics_LCP.i
%include Numerics_AVI.i

%inline %{

  static MixedLinearComplementarityProblem* mixedLinearComplementarityProblemFromFile
    (const char * filename)
  {
    FILE * finput = fopen(filename, "r");
    if (finput)
    {
      MixedLinearComplementarityProblem* problem =
        (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));
      if (mixedLinearComplementarity_newFromFile(problem,finput))
      {
      fprintf(stderr, "mixedLinearComplementarityProblemFromFile: cannot load %s\n",filename);
      free(problem);
      fclose(finput);
      return 0;
      }
      else
      {
        fclose(finput);
        return problem;
      }
    }
    else
    {
      fclose(finput);
      fprintf(stderr, "mixedLinearComplementarityProblemFromFile: cannot open %s\n",filename);
      return 0;
    }
  }
%}

%include Numerics_for_python_callback.i

%include numerics_common.i

%include numerics_MCP.i
%include Numerics_MCP2.i
%include Numerics_VI.i



// FrictionContact
%ignore LocalNonsmoothNewtonSolver; //signature problem (should be SolverOption
                          //instead of *iparam, *dparam).
%ignore DampedLocalNonsmoothNewtonSolver; // signature problem idem.

%ignore frictionContactProblem_new; // signature issue with mu param

%include "typemaps.i"

%apply double *OUTPUT { double *error };
%apply double *OUTPUT { double *result };

// Callback (see SolverOptions.i) needed here
%typemap(in, numinputs=0) (FischerBurmeisterFun3x3Ptr computeACFun3x3) () {
  // Callback (see SolverOptions.i) needed here
  $1 = &fc3d_FischerBurmeisterFunctionGenerated;
 }

%typemap(in, numinputs=0) (AlartCurnierFun3x3Ptr computeACFun3x3) () {
  // Callback (see SolverOptions.i) needed here
  $1 = &fc3d_AlartCurnierFunctionGenerated;
 }

%typemap(in, numinputs=0) (NaturalMapFun3x3Ptr computeACFun3x3) () {
  // Callback (see SolverOptions.i) needed here
  $1 = &fc3d_NaturalMapFunctionGenerated;
 }

// the order matters
%include numerics_FC.i
%include GAMSlink.h
%include numerics_GFC.i

%extend SN_GAMSparams
{

  SN_GAMSparams(SolverOptions* SO)
  {
    assert(SO->solverParameters);
    return (SN_GAMSparams*) SO->solverParameters;
  }

  ~SN_GAMSparams()
  {
    //do nothing
  }

};

//GenericMechanical
//%include GMPReduced.h
//%include GenericMechanicalProblem.h
//%include GenericMechanical_Solvers.h
%{
#include <GenericMechanical_cst.h>
%}
//%include GenericMechanical_cst.h

%include numerics_NM.i

#ifdef WITH_SERIALIZATION
%make_picklable(Callback, Numerics);
%make_picklable(SolverOptions, Numerics);
%make_picklable(FrictionContactProblem, Numerics);
%make_picklable(NumericsMatrix, Numerics);
%make_picklable(SparseBlockStructuredMatrix, Numerics);
#endif
