// NCP
%{
#include "NonlinearComplementarityProblem.h"
#include "NCP_Solvers.h"
#include "NCP_cst.h"
  %}

%include "NonlinearComplementarityProblem.h"
%include "NCP_Solvers.h"
%include "MCP_cst.h"

%extend NonlinearComplementarityProblem
{
  CALL_COMPUTE_F(ncp, (void*))

  CALL_COMPUTE_NABLA_F(ncp, (void*))

  NonlinearComplementarityProblem()
   {
     NonlinearComplementarityProblem* NCP = newNCP();
     return NCP;
   }

  NonlinearComplementarityProblem(SN_OBJ_TYPE* n)
  {
     NonlinearComplementarityProblem* NCP = newNCP();

     SWIG_AsVal_unsigned_SS_int(n, &NCP->n);

     if (NCP->n<1)
     {
       SWIG_Error(SWIG_RuntimeError, "sizeInequalities has to be positive");
       free(NCP);
       return NULL;
     }
     else
     {
     //TODO implement different types of matrices
       NCP->nabla_F = NM_create(NM_DENSE, NCP->n, NCP->n);
     }
     return NCP;
  }


#ifdef SWIGPYTHON
  NonlinearComplementarityProblem(SN_OBJ_TYPE* n, SN_OBJ_TYPE* py_compute)
  {
     NonlinearComplementarityProblem* NCP = newNCP();

     NCP->compute_F = &call_py_compute_Fncp;
     NCP->compute_nabla_F = &call_py_compute_nabla_Fncp;
     SWIG_AsVal_unsigned_SS_int(n, &NCP->n);

     if (NCP->n<1)
     {
       SWIG_Error(SWIG_RuntimeError, "sizeInequalities has to be positive");
       free(NCP);
       return NULL;
     }
     else
     {
       NCP->nabla_F = NM_create(NM_DENSE, NCP->n, NCP->n);
     }

     SN_OBJ_TYPE* method_compute_F = NULL;
     if (PyObject_HasAttrString(py_compute, "compute_F")) method_compute_F = PyObject_GetAttrString(py_compute, "compute_F");
     SN_OBJ_TYPE* method_compute_nabla_F = NULL;
     if (PyObject_HasAttrString(py_compute, "compute_nabla_F")) method_compute_nabla_F = PyObject_GetAttrString(py_compute, "compute_nabla_F");

     if (PyCallable_Check(method_compute_F) && PyCallable_Check(method_compute_nabla_F))
     {
       NCP->env = (void*) malloc(sizeof(class_env_python));
       class_env_python* ncp_env_python = (class_env_python*) NCP->env;
       ncp_env_python->id = ENV_IS_PYTHON_CLASS;
       ncp_env_python->class_object = py_compute;
       target_mem_mgmt_instr(method_compute_F);
       target_mem_mgmt_instr(method_compute_nabla_F);
     }
     else
     {
       target_mem_mgmtX_instr(method_compute_F);
       target_mem_mgmtX_instr(method_compute_nabla_F);
       SWIG_Error(SWIG_TypeError, "argument 2 must be have a method compute_F and a method compute_nabla_F");
       NM_clear(NCP->nabla_F);
       free(NCP->nabla_F);
       free(NCP);
       return NULL;
     }

     return NCP;
   }
#endif /* SWIGPYTHON */

  NonlinearComplementarityProblem(SN_OBJ_TYPE* n, SN_OBJ_TYPE* compute_F, SN_OBJ_TYPE* compute_nabla_F)
  {
     NonlinearComplementarityProblem* NCP = newNCP();

     SWIG_AsVal_unsigned_SS_int(n, &NCP->n);

     if (NCP->n<1)
     {
       SWIG_Error(SWIG_RuntimeError, "sizeInequalities has to be positive");
       free(NCP);
       return NULL;
     }
     else
     {
       NCP->nabla_F = NM_create(NM_DENSE, NCP->n, NCP->n);
     }

     check_save_target_fn(compute_F, NCP->env, env_compute_function, NonlinearComplementarityProblem_call_compute_F, NCP->compute_F, 2);
     check_save_target_fn(compute_nabla_F, NCP->env, env_compute_jacobian, NonlinearComplementarityProblem_call_compute_nabla_F, NCP->compute_nabla_F, 3);

     return NCP;
   }

#ifdef SWIGPYTHON
  void set_compute_F_and_nabla_F_as_C_functions(SN_OBJ_TYPE* lib_name, SN_OBJ_TYPE* compute_F_name, SN_OBJ_TYPE* compute_nabla_F_name)
  {
%#if PY_MAJOR_VERSION < 3
    if(PyString_Check(lib_name) && PyString_Check(compute_F_name) && PyString_Check(compute_nabla_F_name))
%#else
    if(PyUnicode_Check(lib_name) && PyUnicode_Check(compute_F_name) && PyUnicode_Check(compute_nabla_F_name))
%#endif
    {
    void* p_compute_F;
    void* p_compute_nabla_F;

    // TODO: save this lib_handle somewhere and close it !
    get_c_functions(lib_name, compute_F_name, compute_nabla_F_name, &p_compute_F, &p_compute_nabla_F);

    $self->compute_F = (ptrFunctionNCP)p_compute_F;
    $self->compute_nabla_F = (ptrFunctionJacNCP)p_compute_nabla_F;

    }
    else
    {
      SWIG_Error(SWIG_TypeError, "All arguments should be strings");
    }
  }
#endif /* SWIGPYTHON */

    SN_OBJ_TYPE* get_env_as_long(void)
    {
      return SWIG_From_long((uintptr_t)&$self->env);
    }

  ~NonlinearComplementarityProblem()
  {
    if ($self->nabla_F)
    {
      NM_clear($self->nabla_F);
      free($self->nabla_F);
    }
    if ($self->env)
    {
      if(((env_target_lang*)$self->env)->id > 0)
      {
        free($self->env);
      }
    }
    free($self);
  }
};


