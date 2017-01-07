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

  VariationalInequality(PyObject* n)
  {
     VariationalInequality* vi = variationalInequality_new((int) PyInt_AsLong(n));
     vi->F = &call_py_compute_Fvi;
     vi->compute_nabla_F = &call_py_compute_nabla_Fvi;

     if (vi->size < 1)
     {
       PyErr_SetString(PyExc_RuntimeError, "the size of the VI has to be positive");
       free(vi);
       PyErr_PrintEx(0);
       return NULL;
     }

     return vi;
  }


  VariationalInequality(PyObject* n, PyObject* py_compute)
  {

     VariationalInequality* vi = variationalInequality_new((int) PyInt_AsLong(n));

     PyObject* method_compute_F = PyObject_GetAttrString(py_compute, "compute_F");
     PyObject* method_compute_nabla_F = PyObject_GetAttrString(py_compute, "compute_nabla_F");

     if (method_compute_F && method_compute_nabla_F && PyCallable_Check(method_compute_F) && PyCallable_Check(method_compute_nabla_F))
     {
       vi->env = (void*) malloc(sizeof(class_env_python));
       class_env_python* vi_env_python = (class_env_python*) vi->env;
       vi_env_python->id = ENV_IS_PYTHON_CLASS;
       vi_env_python->class_object = py_compute;
       Py_DECREF(method_compute_F);
       Py_DECREF(method_compute_nabla_F);
     }
     else
     {
       if (PyCallable_Check(py_compute))
       {
         vi->F = &call_py_compute_Fvi;
         vi->env = (void*) malloc(sizeof(functions_env_python));
         functions_env_python* vi_env_python = (functions_env_python*) vi->env;
         vi_env_python->id = ENV_IS_PYTHON_FUNCTIONS;
         vi_env_python->env_compute_function = py_compute;
       }
       else
       {
         Py_XDECREF(method_compute_F);
         Py_XDECREF(method_compute_nabla_F);
         PyErr_SetString(PyExc_TypeError, "argument 2 must either be an object with a method compute_F and a method compute_nabla_F or a callable function");
         free(vi);
         PyErr_PrintEx(0);
         return NULL;
       }
     }

     return vi;
   }


   void set_compute_nabla_F(PyObject* py_compute_nabla_F)
   {
     if (PyCallable_Check(py_compute_nabla_F))
     {
       $self->compute_nabla_F = &call_py_compute_nabla_Fvi;
       functions_env_python* vi_env_python = (functions_env_python*) $self->env;
       vi_env_python->id = ENV_IS_PYTHON_FUNCTIONS;
       vi_env_python->env_compute_jacobian = py_compute_nabla_F;
       $self->nabla_F = NM_create(NM_DENSE, $self->size, $self->size);
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 1 must be callable");
       PyErr_PrintEx(0);
     }
   }

   void set_box_constraints(PyObject* box_lower_bound, PyObject* box_upper_bound)
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
         PyErr_PrintEx(0);
         exit(1);
       }
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "The arguments do not have the right length");
       PyErr_PrintEx(0);
       exit(1);
     }
   }

  void set_compute_F_and_nabla_F_as_C_functions(PyObject* lib_name, PyObject* compute_F_name, PyObject* compute_nabla_F_name)
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

    $self->nabla_F = NM_create(NM_DENSE, $self->size, $self->size);
    }
    else
    {
      PyErr_SetString(PyExc_TypeError, "All arguments should be strings");
      PyErr_PrintEx(0);
    }
  }

    PyObject* get_env_as_long(void)
    {
      return PyInt_FromLong((uintptr_t)&$self->env);
    }

  ~VariationalInequality()
  {
    if ($self->nabla_F)
    {
      freeNumericsMatrix($self->nabla_F);
      free($self->nabla_F);
    }
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
          PyErr_SetString(PyExc_TypeError, "unknown set type");
          PyErr_PrintEx(0);
        }
      }
      free($self->set);
      $self->set = NULL;
    }
    if ($self->env)
    {
      if(((env_python*)$self->env)->id > 0)
      {
        free($self->env);
      }
    }
    free($self);
  }
};


