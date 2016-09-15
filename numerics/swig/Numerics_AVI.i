
%include "AffineVariationalInequalities.h"
%include "AVI_Solvers.h"
%include "AVI_cst.h"

%extend AffineVariationalInequalities_
{
  AffineVariationalInequalities_()
   {
     AffineVariationalInequalities_* avi;
     avi = (AffineVariationalInequalities_ *) malloc(sizeof(AffineVariationalInequalities_));
     avi->size = 0;
     avi->M = NULL;
     avi->q = NULL;
     avi->d = NULL;
     avi->poly = NULL;

     return avi;
   }


  AffineVariationalInequalities_(PyObject* mat, PyObject* vec)
  {


     AffineVariationalInequalities_* avi;
     avi = (AffineVariationalInequalities_ *) malloc(sizeof(AffineVariationalInequalities_));
     avi->d = NULL;
     avi->poly = NULL;

      int is_new_object1=0;
      int is_new_object2=0;
      PyArrayObject* array = obj_to_array_fortran_allow_conversion(mat, NPY_DOUBLE,&is_new_object1);
      assert(array);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(vec, NPY_DOUBLE, &is_new_object2); 
      assert(vector);
      NumericsMatrix *M = newNumericsMatrix();
      fillNumericsMatrix(M, NM_DENSE, array_size(array, 0), array_size(array, 1), array_data(array));

      avi->size = M->size0;
      avi->M = M;
      avi->q = (double*)array_data(vector);

      // python mem management
      if(is_new_object1 && array)
      {
        Py_DECREF(array);
      }

      if(is_new_object2 && vector)
      {
        Py_DECREF(vector);
       }

     return avi;
   }

   void set_polyhedron(PyObject* H_mat, PyObject* K_vec)
   {
     $self->poly = (polyhedron*) malloc(sizeof(polyhedron));
      int is_new_object1=0;
      int is_new_object2=0;
      PyArrayObject* array = obj_to_array_fortran_allow_conversion(H_mat, NPY_DOUBLE,&is_new_object1);
      assert(array);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(K_vec, NPY_DOUBLE, &is_new_object2); 
      assert(vector);
      $self->poly->size_ineq = array_size(array, 0);
      $self->poly->size_eq = 0;
      $self->poly->H = (double*)array_data(array);
      $self->poly->K = (double*)array_data(vector);
      $self->poly->Heq = NULL;
      $self->poly->Keq = NULL;

    if (array_size(array, 1) != $self->size)
     {
       PyErr_SetString(PyExc_TypeError, "The matrix does not have the right number of column");
       PyErr_PrintEx(0);
       exit(1);
     }
   }

  ~AffineVariationalInequalities_()
  {
    if ($self->poly)
    {
      free($self->poly);
    }
    if ($self->M)
    {
      free($self->M);
    }
    free($self);
  }
};


