
%include "VariationalInequality.h"
%include "VariationalInequality_Solvers.h"
%include "VI_cst.h"

%extend VariationalInequality
{
  VariationalInequality()
   {
     VariationalInequality* vi;
     vi = (VariationalInequality *) malloc(sizeof(VariationalInequality));
     vi->nabla_F = NULL;
     vi->F = &call_py_compute_F;
     vi->compute_nabla_F = &call_py_compute_nabla_F;
     vi->env = NULL;

     return vi;
   }

  VariationalInequality(PyObject* n, PyObject* py_compute)
  {
     VariationalInequality* vi;
     vi =  (VariationalInequality *) malloc(sizeof(VariationalInequality));

     vi->size = (int) PyInt_AsLong(n);
     vi->set = NULL;
     vi->env = NULL;
     vi->compute_nabla_F = NULL;
     vi->ProjectionOnX = NULL;
     vi->nabla_F = NULL;
     vi->normVI = std::numeric_limits<double>::quiet_NaN();
     vi->istheNormVIset = 0;


     if (vi->size < 1)
     {
       PyErr_SetString(PyExc_RuntimeError, "the size of the VI has to be positive");
       free(vi);
       return NULL;
     }


     PyObject* method_compute_F = PyObject_GetAttrString(py_compute, "compute_F");
     PyObject* method_compute_nabla_F = PyObject_GetAttrString(py_compute, "compute_nabla_F");

     if (PyCallable_Check(method_compute_F) && PyCallable_Check(method_compute_nabla_F))
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
         vi->F = &call_py_compute_F;
         vi->env = (void*) malloc(sizeof(functions_env_python));
         functions_env_python* vi_env_python = (functions_env_python*) vi->env;
         vi_env_python->id = -1;
         vi_env_python->env_compute_function = py_compute;
       }
       else
       {
         Py_XDECREF(method_compute_F);
         Py_XDECREF(method_compute_nabla_F);
         PyErr_SetString(PyExc_TypeError, "argument 2 must either be an object with a method compute_F and a method compute_nabla_F or a callable function");
         free(vi->nabla_F);
         free(vi);
         PyErr_PrintEx(0);
         exit(1);
       }
     }

     return vi;
   }


   void set_compute_nabla_F(PyObject* py_compute_nabla_F)
   {
     if (PyCallable_Check(py_compute_nabla_F))
     {
       $self->compute_nabla_F = &call_py_compute_nabla_F;
       functions_env_python* vi_env_python = (functions_env_python*) $self->env;
       vi_env_python->id = ENV_IS_PYTHON_FUNCTIONS;
       vi_env_python->env_compute_jacobian = py_compute_nabla_F;
       $self->nabla_F = (double*) calloc($self->size*$self->size, sizeof(double));
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
       box_c->id = VI_SET_IS_BOX;
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



  ~VariationalInequality()
  {
    if ($self->nabla_F)
    {
      free($self->nabla_F);
    }
    if ($self->set)
    {
      //black magic
      switch (((generic_set*)$self->set)->id)
      {
        case VI_SET_IS_BOX:
        {
          box_constraints* box_c = (box_constraints*)$self->set;
          free_box(box_c);
          break;
        }
        default:
          PyErr_SetString(PyExc_TypeError, "unknown set type");
      }
      free($self->set);
      $self->set = NULL;
    }
    free($self);
  }
};


