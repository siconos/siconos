
%apply (double *b) { (double *b_bck_gfc3d) };

// for b in GlobalFrictionContact
%typemap(out) (double* b) {
  npy_intp dims[1];

  if (!arg1->H) { PyErr_SetString(PyExc_TypeError, "H is not present, don't known the size"); SWIG_fail; }

  dims[0] = arg1->H->size1;
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

// vectors of size col(H)
%typemap(memberin) (double *b) {
  // Still some dark magic :( --xhub
 char msg[1024];
  assert(arg1);
  if (!arg1->H)
  {
    PyErr_SetString(PyExc_RuntimeError, "H is not initialized, it sould be done first!");
    SWIG_fail;
  }

  int size = arg1->H->size1;
  if (size !=  array_size(array2, 0))
  {
    snprintf(msg, sizeof(msg), "Size of b is %ld, but the size of H is %d! Both should be equal!\n", array_size(array2, 0), size);
    PyErr_SetString(PyExc_RuntimeError, msg);
    SWIG_fail;
  }

  if (!$1) { $1 = (double*)malloc(size * sizeof(double)); }
  memcpy($1, $input, size * sizeof(double));

}

%{
#include "GlobalFrictionContactProblem.h"
#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Solvers.h"
#include "gfc3d_compute_error.h"
  
%}

%include "GlobalFrictionContactProblem.h"
%include "gfc3d_Solvers.h"
%include "gfc3d_compute_error.h"

%extend GlobalFrictionContactProblem
{

  GlobalFrictionContactProblem()
    {

      GlobalFrictionContactProblem *FC;
      // return pointer : free by std swig destructor
      FC = (GlobalFrictionContactProblem *) malloc(sizeof(GlobalFrictionContactProblem));

      globalFrictionContact_null(FC);

      return FC;
    }

  GlobalFrictionContactProblem(PyObject *dim)
    {

      GlobalFrictionContactProblem *FC;
      // return pointer : free by std swig destructor
      FC = (GlobalFrictionContactProblem *) malloc(sizeof(GlobalFrictionContactProblem));

      globalFrictionContact_null(FC);
      FC->dimension = (int) PyInt_AsLong(dim);

      return FC;
    }

  GlobalFrictionContactProblem(PyObject *dim, PyObject *o1, PyObject *o2, PyObject *o3)
    {

      int is_new_object1=0;
      int is_new_object2=0;
      int is_new_object3=0;

      PyArrayObject* array = obj_to_array_fortran_allow_conversion(o1, NPY_DOUBLE,&is_new_object1);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(o2, NPY_DOUBLE, &is_new_object2);
      PyArrayObject* mu_vector = obj_to_array_contiguous_allow_conversion(o3, NPY_DOUBLE, &is_new_object3);
      GlobalFrictionContactProblem * FC = (GlobalFrictionContactProblem *) malloc(sizeof(GlobalFrictionContactProblem));
      globalFrictionContact_null(FC);

      size_t size0 = array_size(array,0);
      size_t size1 = array_size(array,1);
      FC->M = NM_create(NM_DENSE, size0, size1);

      memcpy(FC->M->matrix0,array_data(array),size0*size1*sizeof(double));
      FC->dimension = (int) PyInt_AsLong(dim);
      FC->numberOfContacts = size0 / FC->dimension;
      FC->q = (double *) malloc(size0*sizeof(double));
      memcpy(FC->q,array_data(vector),size0*sizeof(double));
      FC->mu = (double *) malloc(FC->numberOfContacts*sizeof(double));
      memcpy(FC->mu,array_data(mu_vector),FC->numberOfContacts*sizeof(double));

      // python mem management
      if(is_new_object1 && array)
      {
        Py_DECREF(array);
      }

      if(is_new_object2 && vector)
      {
        Py_DECREF(vector);
      }

      if(is_new_object3 && mu_vector)
      {
        Py_DECREF(mu_vector);
      }

      return FC;
    }

  ~GlobalFrictionContactProblem()
  {
    freeGlobalFrictionContactProblem($self);
  }

};


// remove typemaps for b
%clear double* b;
%apply (double *b_bck_gfc3d) { (double *b) };

