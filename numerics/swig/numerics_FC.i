%{
#include "FrictionContactProblem.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "Friction_cst.h"
#ifdef WITH_FCLIB
#include "fclib_interface.h"
#endif
%}

%include "FrictionContactProblem.h"
#ifdef WITH_FCLIB
%include fclib_interface.h
#endif

%include "fc3d_Solvers.h"
%include "fc3d_unitary_enumerative.h"
%include "fc2d_Solvers.h"
%include "Friction_cst.h"
%include "fc3d_compute_error.h"

%extend FrictionContactProblem_
{
  FrictionContactProblem_()
  {
    FrictionContactProblem * FCP = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    FCP->M = NULL;
    FCP->q = NULL;
    FCP->mu = NULL;

    return FCP;
  }


  /* copy constructor */
  FrictionContactProblem_(PyObject *o)
  {
    FrictionContactProblem* fcp;
    FrictionContactProblem* FCP;

    %SN_INPUT_CHECK_RETURN(o, fcp, FrictionContactProblem);

    FCP = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    FCP->dimension = fcp->dimension;
    FCP->M = fcp->M;
    FCP->numberOfContacts = fcp->numberOfContacts;
    FCP->q =  fcp->q;
    FCP->mu = fcp->mu;

    Py_INCREF(o);

    return FCP;
  }

  /* */
  FrictionContactProblem_(PyObject *dim, PyObject *numberOfContacts, PyObject *M, PyObject *q, PyObject *mu)
  {
    FrictionContactProblem * FC = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
    FC->dimension = PyInt_AsLong(dim);
    FC->numberOfContacts = PyInt_AsLong(numberOfContacts);

    %SN_INPUT_CHECK_RETURN(q, FC->M, NumericsMatrix);
    %SN_INPUT_CHECK_RETURN(q, FC->q, double);
    %SN_INPUT_CHECK_RETURN(mu, FC->mu, double);

    return FC;

  }

  NumericsMatrix* rawM()
  {
    return $self->M;
  }

  FrictionContactProblem_(PyObject *dim, PyObject *o1, PyObject *o2, PyObject *o3)
    {

      int is_new_object1=0;
      int is_new_object2=0;
      int is_new_object3=0;

      PyArrayObject* array = obj_to_array_fortran_allow_conversion(o1, NPY_DOUBLE,&is_new_object1);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(o2, NPY_DOUBLE, &is_new_object2);
      PyArrayObject* mu_vector = obj_to_array_contiguous_allow_conversion(o3, NPY_DOUBLE, &is_new_object3);
      FrictionContactProblem *FC;
      // return pointer : free by std swig destructor
      FC = (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
      size_t size0 = array_size(array,0);
      size_t size1 = array_size(array,1);

      FC->M = createNumericsMatrix(NM_DENSE, size0, size1);

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

  ~FrictionContactProblem_()
  {
    freeFrictionContactProblem($self);
  }

};

%inline %{

#include <stdio.h>
  static FrictionContactProblem* frictionContactProblemFromFile
    (const char * filename)
  {
    FILE * finput = fopen(filename, "r");
    if (finput)
    {
      FrictionContactProblem* problem =
        (FrictionContactProblem *) malloc(sizeof(FrictionContactProblem));
      if (frictionContact_newFromFile(problem,finput))
      {
      char msg[1024];
      snprintf(msg, sizeof(msg), "frictionContactProblemFromFile: cannot load %s\n",filename);
      PyErr_SetString(PyExc_RuntimeError, msg);
      PyErr_PrintEx(0);
      free(problem);
      fclose(finput);
      return NULL;
      }
      else
      {
        fclose(finput);
        return problem;
      }
    }
    else
    {
      char msg[1024];
      snprintf(msg, sizeof(msg), "frictionContactProblemFromFile: cannot open %s\n",filename);
      PyErr_SetString(PyExc_RuntimeError, msg);
      PyErr_PrintEx(0);
      return NULL;
    }
    
  }

%}
