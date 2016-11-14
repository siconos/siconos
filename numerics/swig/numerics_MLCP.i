// MLCP
%apply (double *q) { (double *qbck) };

// redefine typemap on q for MLCP
%typemap(out) (double* q) {
  npy_intp dims[1];

  dims[0] = arg1->m + arg1->n;
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

%{
  
  #include "mlcp_cst.h"
  #include "MixedLinearComplementarityProblem.h"
  #include "MLCP_Solvers.h"
  
  %}

%include "mlcp_cst.h"
%include "MixedLinearComplementarityProblem.h"
%include "MLCP_Solvers.h"

%exception MixedLinearComplementarityProblem {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%extend MixedLinearComplementarityProblem
{
  MixedLinearComplementarityProblem()
   {
     MixedLinearComplementarityProblem* MLCP;
     MLCP =  (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));
     return MLCP;
   }

  MixedLinearComplementarityProblem(PyObject *dim, PyObject *o1, PyObject *o2)
    {

      int is_new_object1=0;
      int is_new_object2=0;

      PyArrayObject* array = obj_to_array_fortran_allow_conversion(o1, NPY_DOUBLE,&is_new_object1);
      PyArrayObject* vector = obj_to_array_contiguous_allow_conversion(o2, NPY_DOUBLE, &is_new_object2);

      if ( array_size(array,0) !=  array_size(array,1))
      {
        PyErr_Format(PyExc_ValueError,
                     "A non square matrix (%ld,%ld) has been given",
                     array_size(array,0), array_size(array,1));
      }


      MixedLinearComplementarityProblem *MLCP;
      // return pointer : free by std swig destructor
      MLCP = (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));

      MLCP->M = NM_create(NM_DENSE, array_size(array,0), array_size(array,1));

      memcpy(MLCP->M->matrix0,array_data(array),MLCP->M->size0*MLCP->M->size1*sizeof(double));

      MLCP->n = (int) PyInt_AsLong(dim);
      MLCP->m = MLCP->M->size0 - MLCP->n;
      MLCP->blocksRows = (int *) malloc(3*sizeof(int));
      MLCP->blocksIsComp = (int *) malloc(2*sizeof(int));


      MLCP->blocksRows[0]=0;
      MLCP->blocksRows[1]=MLCP->n;
      MLCP->blocksRows[2]=MLCP->n+MLCP->m;
      MLCP->blocksIsComp[0]=0;
      MLCP->blocksIsComp[1]=1;


      MLCP->isStorageType1 = 1;
      MLCP->isStorageType2 = 0;
      MLCP->A = NULL;
      MLCP->B = NULL;
      MLCP->C = NULL;
      MLCP->D = NULL;
      MLCP->a = NULL;
      MLCP->b = NULL;

      if ( array_size(array,0) !=  array_size(vector,0))
      {
        //printf("size of q = %i\n",  array_size(vector,0));
        //printf("size of M = %i\n",  array_size(array,0));

        PyErr_Format(PyExc_ValueError,
                     "Matrix and vector of incompatible lengths (%ld != %ld) ",
                     array_size(array,0), array_size(vector,0) );
      }
      MLCP->q = (double *) malloc(MLCP->M->size0*sizeof(double));
      memcpy(MLCP->q,array_data(vector),MLCP->M->size0*sizeof(double));

      // python mem management
      if(is_new_object1 && array)
      {
        Py_DECREF(array);
      }

      if(is_new_object2 && vector)
      {
        Py_DECREF(vector);
      }

      return MLCP;

    }



  ~MixedLinearComplementarityProblem()
  {
    freeMixedLinearComplementarityProblem($self);
  }

  // MixedLinearComplementarityProblem * newFromFilename(PyObject * o1)
  // {
  //   int result;
  //   MixedLinearComplementarityProblem *MLCP;
  //   // return pointer : free by std swig destructor
  //   MLCP = (MixedLinearComplementarityProblem *) malloc(sizeof(MixedLinearComplementarityProblem));

  //   char *arg1 = (char *) 0 ;
  //   int res1 ;
  //   char *buf1 = 0 ;
  //   int alloc1 = 0 ;

  //   res1 = SWIG_AsCharPtrAndSize(o1, &buf1, NULL, &alloc1);
  //   // if (!SWIG_IsOK(res1)) {
  //   //   SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "MixedLinearComplementarity_newFromFilename" "', argument " "1"" of type '" "char *""'");
  //   // }
  //   arg1 = reinterpret_cast< char * >(buf1);
  //   {
  //     try
  //     {
  //       result = (int)mixedLinearComplementarity_newFromFilename(MLCP,arg1);
  //     }
  //     catch (const std::invalid_argument& e)
  //     {
  //       // SWIG_exception(SWIG_ValueError, e.what());
  //     }
  //   }

  //   return MLCP;

  // }

};


%clear double* q;
%apply (double *qbck) { (double *q) };
