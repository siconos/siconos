// LCP

%ignore lcp_compute_error_only;

%{
  #include "LinearComplementarityProblem.h"
  #include "LCP_Solvers.h"
  #include "lcp_cst.h"
  %}

%include "LinearComplementarityProblem.h"
%include "LCP_Solvers.h"
%include "lcp_cst.h"

%extend LinearComplementarityProblem
{
  LinearComplementarityProblem(PyObject *o1, PyObject *o2)
    {

      int is_new_object1=0;
      int is_new_object2=0;
      PyArrayObject* array = obj_to_array_fortran_allow_conversion(o1, NPY_DOUBLE,&is_new_object1);
      assert(array);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(o2, NPY_DOUBLE, &is_new_object2);
      assert(vector);
      LinearComplementarityProblem *LC;
      // return pointer : free by std swig destructor
      LC = (LinearComplementarityProblem *) malloc(sizeof(LinearComplementarityProblem));
      size_t size0 = array_size(array,0);
      size_t size1 = array_size(array,1);
      LC->M = NM_create(NM_DENSE, size0, size1);

      memcpy(LC->M->matrix0,array_data(array),size0*size1*sizeof(double));
      LC->size = size0;
      LC->q = (double *) malloc(size0*sizeof(double));
      memcpy(LC->q,array_data(vector),size0*sizeof(double));

      // python mem management
      if(is_new_object1 && array)
      {
        Py_DECREF(array);
      }

      if(is_new_object2 && vector)
      {
        Py_DECREF(vector);
       }

      return LC;

    }


  ~LinearComplementarityProblem()
  {
    freeLinearComplementarityProblem($self);
  }

};



