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
  NonlinearComplementarityProblem()
   {
     NonlinearComplementarityProblem* NCP;
     NCP = (NonlinearComplementarityProblem *) malloc(sizeof(NonlinearComplementarityProblem));
     NCP->nabla_F = NULL;
     NCP->compute_F = &call_py_compute_Fncp;
     NCP->compute_nabla_F = &call_py_compute_nabla_Fncp;
     NCP->env = NULL;

     return NCP;
   }

  NonlinearComplementarityProblem(PyObject* n)
  {
     NonlinearComplementarityProblem* NCP;
     NCP =  (NonlinearComplementarityProblem *) malloc(sizeof(NonlinearComplementarityProblem));

     NCP->compute_F = &call_py_compute_Fncp;
     NCP->compute_nabla_F = &call_py_compute_nabla_Fncp;
     NCP->n = (int) PyInt_AsLong(n);

     if (NCP->n<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeInequalities has to be positive");
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


  NonlinearComplementarityProblem(PyObject* n, PyObject* py_compute_class)
  {
     NonlinearComplementarityProblem* NCP;
     NCP =  (NonlinearComplementarityProblem *) malloc(sizeof(NonlinearComplementarityProblem));

     NCP->compute_F = &call_py_compute_Fncp;
     NCP->compute_nabla_F = &call_py_compute_nabla_Fncp;
     NCP->n = (int) PyInt_AsLong(n);

     if (NCP->n<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeInequalities has to be positive");
       free(NCP);
       return NULL;
     }
     else
     {
       NCP->nabla_F = NM_create(NM_DENSE, NCP->n, NCP->n);
     }

     PyObject* method_compute_F = PyObject_GetAttrString(py_compute_class, "compute_F");
     PyObject* method_compute_nabla_F = PyObject_GetAttrString(py_compute_class, "compute_nabla_F");
     if (PyCallable_Check(method_compute_F) && PyCallable_Check(method_compute_nabla_F))
     {
       NCP->env = (void*) malloc(sizeof(class_env_python));
       class_env_python* ncp_env_python = (class_env_python*) NCP->env;
       ncp_env_python->id = ENV_IS_PYTHON_CLASS;
       ncp_env_python->class_object = py_compute_class;
       Py_DECREF(method_compute_F);
       Py_DECREF(method_compute_nabla_F);
     }
     else
     {
       Py_XDECREF(method_compute_F);
       Py_XDECREF(method_compute_nabla_F);
       PyErr_SetString(PyExc_TypeError, "argument 2 must be have a method compute_F and a method compute_nabla_F");
       freeNumericsMatrix(NCP->nabla_F);
       free(NCP->nabla_F);
       free(NCP);
       PyErr_PrintEx(0);
       exit(1);
     }

     return NCP;
   }

  NonlinearComplementarityProblem(PyObject* n, PyObject* py_compute_F, PyObject* py_compute_nabla_F)
  {
     NonlinearComplementarityProblem* NCP;
     NCP =  (NonlinearComplementarityProblem *) malloc(sizeof(NonlinearComplementarityProblem));

     NCP->compute_F = &call_py_compute_Fncp;
     NCP->compute_nabla_F = &call_py_compute_nabla_Fncp;
     NCP->n = (int) PyInt_AsLong(n);

     if (NCP->n<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeInequalities has to be positive");
       free(NCP);
       return NULL;
     }
     else
     {
       NCP->nabla_F = NM_create(NM_DENSE, NCP->n, NCP->n);
     }

     NCP->env = (void*) malloc(sizeof(functions_env_python));
     functions_env_python* ncp_env_python = (functions_env_python*) NCP->env;
     ncp_env_python->id = ENV_IS_PYTHON_FUNCTIONS;

     if (PyCallable_Check(py_compute_F)) 
     {
       ncp_env_python->env_compute_function = py_compute_F;
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 3 must be callable");
       freeNumericsMatrix(NCP->nabla_F);
       free(NCP->nabla_F);
       free(NCP->env);
       free(NCP);
       PyErr_PrintEx(0);
       exit(1);
     }


     if (PyCallable_Check(py_compute_nabla_F))
     {
       ncp_env_python->env_compute_jacobian = py_compute_nabla_F;
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 4 must be callable");
       freeNumericsMatrix(NCP->nabla_F);
       free(NCP->nabla_F);
       free(NCP->env);
       free(NCP);
       PyErr_PrintEx(0);
       exit(1);
     }

     return NCP;
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

    $self->compute_F = (ptrFunctionNCP)p_compute_F;
    $self->compute_nabla_F = (ptrFunctionJacNCP)p_compute_nabla_F;

    }
    else
    {
      PyErr_SetString(PyExc_TypeError, "All arguments should be strings");
    }
  }

    PyObject* get_env_as_long(void)
    {
      return PyInt_FromLong((uintptr_t)&$self->env);
    }

  ~NonlinearComplementarityProblem()
  {
    if ($self->nabla_F)
    {
      freeNumericsMatrix($self->nabla_F);
      free($self->nabla_F);
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


