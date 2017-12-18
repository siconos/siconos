%inline %{
#include "VariationalInequality.h"
#include "VariationalInequality_Solvers.h"
#include "VI_cst.h"
#include "VariationalInequality_Solvers.h"
%}

%include "VariationalInequality.h"
%include "VariationalInequality_Solvers.h"
%include "VI_cst.h"

%extend VariationalInequality
{

  CALL_COMPUTE_F(vi, VI_get_env)

  CALL_COMPUTE_NABLA_F(vi, VI_get_env)

  VariationalInequality(SN_OBJ_TYPE* n)
  {
     int nn;
     SWIG_AsVal_int(n, &nn);
     VariationalInequality* vi = variationalInequality_new(nn);
     assert(vi);

     if (vi->size < 1)
     {
       SWIG_Error(SWIG_RuntimeError, "the size of the VI has to be positive");
       TARGET_ERROR_VERBOSE;
       free(vi);
       return NULL;
     }

     vi->nabla_F = NM_create(NM_DENSE, nn, nn);

     return vi;
  }


  VariationalInequality(SN_OBJ_TYPE* n, SN_OBJ_TYPE* compute)
  {

     int nn;
     SWIG_AsVal_int(n, &nn);
     VariationalInequality* vi = variationalInequality_new(nn);

     vi->nabla_F = NM_create(NM_DENSE, nn, nn);

#ifdef SWIGPYTHON
     SN_OBJ_TYPE* method_compute_F = NULL;
     if (PyObject_HasAttrString(compute, "compute_F")) method_compute_F = PyObject_GetAttrString(compute, "compute_F");
     SN_OBJ_TYPE* method_compute_nabla_F = NULL;
     if (PyObject_HasAttrString(compute, "compute_nabla_F")) method_compute_nabla_F = PyObject_GetAttrString(compute, "compute_nabla_F");

     if (method_compute_F && method_compute_nabla_F && PyCallable_Check(method_compute_F) && PyCallable_Check(method_compute_nabla_F))
     {
       vi->env = (void*) malloc(sizeof(class_env_python));
       class_env_python* vi_env_python = (class_env_python*) vi->env;
       vi_env_python->id = ENV_IS_PYTHON_CLASS;
       vi_env_python->class_object = compute;
       target_mem_mgmt_instr(method_compute_F);
       target_mem_mgmt_instr(method_compute_nabla_F);
     }
     else
#endif /* SWIGPYTHON */
     {
       check_save_target_fn(compute, vi->env, env_compute_function, VariationalInequality_call_compute_F, vi->F, 2);
     }

#ifdef SWIGPYTHON
     target_mem_mgmtX_instr(method_compute_F);
     target_mem_mgmtX_instr(method_compute_nabla_F);
#endif /* SWIGPYTHON */
     return vi;
   }


   void set_compute_nabla_F(SN_OBJ_TYPE* compute_nabla_F)
   {
     check_save_target_fn(compute_nabla_F, $self->env, env_compute_jacobian, VariationalInequality_call_compute_nabla_F, $self->compute_nabla_F, 1);
   }

#ifdef SWIGPYTHON
   void set_box_constraints(SN_OBJ_TYPE* box_lower_bound, SN_OBJ_TYPE* box_upper_bound)
   {
     if ((PyObject_Length(box_lower_bound) == $self->size) && (PyObject_Length(box_upper_bound) == $self->size))
     {
       $self->set = (void *)malloc(sizeof(box_constraints));
       box_constraints* box_c = (box_constraints*)$self->set;
       box_c->id = SICONOS_SET_BOX;
       box_c->lb = (double*)malloc($self->size*sizeof(double));
       box_c->ub = (double*)malloc($self->size*sizeof(double));

       if (!convert_darray(box_lower_bound, box_c->lb) || !convert_darray(box_upper_bound, box_c->ub))
       {
         TARGET_ERROR_VERBOSE;
         exit(1);
       }
     }
     else
     {
       SWIG_Error(SWIG_TypeError, "The arguments do not have the right length");
       TARGET_ERROR_VERBOSE;
     }
   }

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

    $self->F = (ptrFunctionVI)p_compute_F;
    $self->compute_nabla_F = (ptrFunctionVI_nabla)p_compute_nabla_F;
    }
    else
    {
      SWIG_Error(SWIG_TypeError, "All arguments should be strings");
      TARGET_ERROR_VERBOSE;
    }
  }
#endif /* SWIGPYTHON */

    SN_OBJ_TYPE* get_env_as_long(void)
    {
      return SWIG_From_long((uintptr_t)&$self->env);
    }

  ~VariationalInequality()
  {
    if ($self->set)
    {
      //black magic
      switch (((generic_set*)$self->set)->id)
      {
        case SICONOS_SET_BOX:
        {
          box_constraints* box_c = (box_constraints*)$self->set;
          free_box(box_c);
          break;
        }
        default:
        {
          SWIG_Error(SWIG_TypeError, "unknown set type");
          TARGET_ERROR_VERBOSE;
        }
      }
      free($self->set);
      $self->set = NULL;
    }
    if ($self->env)
    {
      if(((env_target_lang*)$self->env)->id > 0)
      {
        free($self->env);
      }
    }
    freeVariationalInequalityProblem($self);
  }
};


