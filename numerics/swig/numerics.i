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

%include exception.i

%{
#include "tlsdef.h"
#include "sn_error_handling.h"

tlsvar char error_msg[2048];

static char* format_exception_msg(const char* first_line)
{
  strncpy(error_msg, first_line, strlen(first_line)+1);
  strncat(error_msg, "\n", 2);
  const char* sn_msg = sn_fatal_error_msg();
  strncat(error_msg, sn_msg, strlen(sn_msg) - 1);
  return error_msg;
}

static char* format_msg_concat(const char* msg1, const char* msg2)
{
  strncpy(error_msg, msg1, strlen(msg1)+1);
  strncat(error_msg, "\n", 2);
  strncat(error_msg, msg2, strlen(msg2));
  return error_msg;
}

%}

%exception {
/* I'm HERE */
/* TODO implement setjmp/longjmp here  + SWIG_exception */
  switch (SN_SETJMP_EXTERNAL_START)
  {
  case SN_NO_ERROR:
  {
    $action
    SN_SETJMP_EXTERNAL_STOP
    break;
  }
  case SN_MEMORY_ALLOC_ERROR:
  {
    SWIG_exception(SWIG_MemoryError, format_exception_msg("Out of memory:"));
    break;
  }
  case SN_UNSUPPORTED_LINALG_OP:
  {
    SWIG_exception(SWIG_RuntimeError, format_exception_msg("Unsupported linear algebra operation:"));
    break;
  }
  case SN_PROBLEM_NOT_PROCESSABLE:
  {
    SWIG_exception(SWIG_RuntimeError, format_exception_msg("The given problem is not processable:"));
    break;
  }
  default:
  {
    SWIG_exception(SWIG_UnknownError, format_exception_msg("Unknown error! Hopefully more info follow:"));
    break;
  }
  }

}

// generated docstrings from doxygen xml output
%include numerics-docstrings.i

// list of functions that returns objects that they have created
%include numerics_newobjects.i

%include <typemaps/swigmacros.swg>

%define %SN_INPUT_CHECK_RETURN(var, output, type_)
{ void* _ptr = NULL; if (!SWIG_IsOK(SWIG_ConvertPtr(var, &_ptr, SWIGTYPE_p_##type_, 0 |  0 ))) {
    char errmsg[1024];
    snprintf(errmsg, sizeof(errmsg), "Argument check failed! Argument %s has the wrong type, should be %s", #var, #type_);
    SWIG_Error(SWIG_TypeError, errmsg); return NULL; }
    output = %static_cast(_ptr, type_*);
   }
%enddef

%{
#include "SiconosConfig.h"
#include "SiconosNumerics.h"
#include "SolverOptions.h"
#include "SparseMatrix_internal.h"
#include "NumericsMatrix.h"
#include "SparseBlockMatrix.h"
#include "NumericsSparseMatrix.h"

#ifdef SWIGPYTHON
#include "Numerics_python_functions.h"
#endif /* SWIGPYTHON */

#ifdef SWIGMATLAB
#include "Numerics_matlab_functions.h"
#endif /* SWIGMATLAB */

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
#ifdef __cplusplus
#if defined(SICONOS_STD_SHARED_PTR) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
namespace std11 = std;
#include <memory>
#else
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
namespace std11 = boost;
#endif
#endif
%}
#ifdef __cplusplus
#define SWIG_SHARED_PTR_NAMESPACE std11
%include boost_shared_ptr.i
#endif

// put needed shared_ptrs here
// commented-out for now, swig insists on calling freeFrictionContactProblem()
// instead of respecting the shared_ptr!
// %shared_ptr(FrictionContactProblem)


 // more convenient
 %rename (LCP) LinearComplementarityProblem;
 %rename (MLCP) MixedLinearComplementarityProblem;
 %rename (MCP) MixedComplementarityProblem;
 %rename (NCP) NonlinearComplementarityProblem;
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
%import tlsdef.h
%include numerics_verbose.h

// this has to die --xhub
// info param
%typemap(in, numinputs=0) (int *info) (int temp_info = -1)
{
  // a default initialization : solver may stop if *info = 0 (checkTrivialCase)
  // checkTrivialCase => better if directly in solvers, not in driver.
  $1 = &temp_info;
}
%warnfilter(322) set_cstruct;
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

#ifdef SWIGPYTHON
%fragment("NumPy_Fragments");
#endif /* SWIGPYTHON */

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
#include "NCP_cst.h"
#include "mlcp_cst.h"
#include "VI_cst.h"
#include "ConvexQP_cst.h"
#include "GenericMechanical_cst.h"
#include "fc2d_Solvers.h"
#include "fc3d_Solvers.h"
#include "gfc3d_Solvers.h"
#include "MCP_Solvers.h"
#include "NCP_Solvers.h"
#include "MLCP_Solvers.h"
#include "ConvexQP_Solvers.h"
#include "SiconosCompat.h"
#include "SOCLCP_Solvers.h"
#include "NonSmoothDrivers.h"
  %}

%include numerics_NM.i

%include numerics_MLCP.i

/////////////////////////
// This is common to all the problem defined below
// MLCP should be defined above this statement
/////////////////////////

%typemap(out) (double* q) {

  if (!arg1->M) { SWIG_exception_fail(SWIG_RuntimeError, "M is not present, don't known the size"); }

  if ($1)
  {
    SN_OBJ_TYPE *obj;
    C_to_target_lang1(obj, arg1->M->size0, $1, SWIG_fail);
    $result = obj;
  }
  else
    SWIG_fail;
 }

%typemap(out) (double* mu) {

  if (arg1->numberOfContacts <= 0) { SWIG_exception_fail(SWIG_RuntimeError, "numberOfContacts is not set"); }

  if ($1)
  {
    SN_OBJ_TYPE *obj;
    C_to_target_lang1(obj, arg1->numberOfContacts, $1, SWIG_fail);
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
    SWIG_exception(SWIG_RuntimeError, "numberOfContacts is not set, it sould be done first!");
    SWIG_fail;
  }

  if (arg1->numberOfContacts !=  array_size(array2, 0))
  {
    char msg[1024];
    snprintf(msg, sizeof(msg), "Size of mu is %ld, but the number of contacts is %d! Both should be equal!\n", array_size(array2, 0), arg1->numberOfContacts);
    SWIG_exception_fail(SWIG_ValueError, msg);
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
    SWIG_exception_fail(SWIG_RuntimeError, "M is not initialized, it sould be done first!");
  }

  int size = arg1->M->size0;
  if (size !=  array_size(array2, 0))
  {
    snprintf(msg, sizeof(msg), "Size of q is %ld, but the size of M is %d! Both should be equal!\n", array_size(array2, 0), size);
    SWIG_exception_fail(SWIG_RuntimeError, msg);
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


%include numerics_common.i

%include Numerics_callback.i

#ifdef SWIGPYTHON
%include Numerics_for_python_callback.i
%include numerics_MCP.i
#endif /* SWIGPYTHON */
%include Numerics_MCP2.i
%include Numerics_NCP.i
%include Numerics_VI.i
%include Numerics_ConvexQP.i
%include Numerics_SOCLCP.i



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

%define STR_FIELD_COPY(field,strobj)
{
  int alloc = 0;
  char* name_str;
  size_t len = 0;
  int res = SWIG_AsCharPtrAndSize(strobj, &name_str, &len, &alloc);
  if (!SWIG_IsOK(res)) {
    SWIG_Error(SWIG_ArgError(res), "in method unknown', argument " "1"" of type '" "char *""'");
  }

  // Some weird logic here
  if (field)
  {
    field = (char*)realloc(field, len*sizeof(char));
  }
  else
  {
    field = (char*)malloc(len*sizeof(char));
  }
  strncpy(field, name_str, len);

  if (alloc == SWIG_NEWOBJ) free(name_str);

}
%enddef

%extend SN_GAMSparams
{

  SN_GAMSparams(SolverOptions* SO)
  {
    assert(SO->solverParameters);
    return (SN_GAMSparams*) SO->solverParameters;
  }

  void gamsdir_set(SN_OBJ_TYPE* strobj)
  {
    STR_FIELD_COPY($self->gams_dir, strobj)
  }

  void modeldir_set(SN_OBJ_TYPE* strobj)
  {
    STR_FIELD_COPY($self->model_dir, strobj)
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

#ifdef WITH_SERIALIZATION
%make_picklable(Callback, Numerics);
%make_picklable(SolverOptions, Numerics);
%make_picklable(FrictionContactProblem, Numerics);
%make_picklable(NumericsMatrix, Numerics);
%make_picklable(SparseBlockStructuredMatrix, Numerics);
#endif
