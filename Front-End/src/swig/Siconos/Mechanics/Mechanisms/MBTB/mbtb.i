/* mbtb.i this file contains exported API of the MBTB library.*/
%module mbtb
%{
#define SWIG_FILE_WITH_INIT
#include <sstream>
#if defined(Py_COMPLEXOBJECT_H)
#undef c_sum
#undef c_diff
#undef c_neg
#undef c_prod
#undef c_quot
#undef c_pow
#undef c_abs
#endif
%}
%{
#include "MBTB_PYTHON_API.hpp"
  //#include "MBTB_Body.hpp"
%}


// mandatory !
%rename (lambda_) lambda;

// shared ptr management
%include "boost_shared_ptr.i"

// numpy macros
%include numpy.i 	

%init %{
  import_array();
%}

// handle standard exceptions
%include "exception.i"
%exception
{
  try
  {
    $action
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch (const std::out_of_range& e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
}

// handle stl data types
%include "stl.i"

// 1. Vector and Matrix <=> numpy array (dense only)

%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosVector>)
{
  int res = SWIG_ConvertPtr($input, 0, SWIGTYPE_p_boost__shared_ptrT_SiconosVector_t, 0);
  int _v = SWIG_CheckState(res);
  $1 = is_array($input) || PySequence_Check($input) || _v;
}
%typemap(in,fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> (PyArrayObject* array=NULL, int is_new_object)
{

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array)
  {
    void *argp;
    SWIG_fail; // not implemented : $1 = type_conv($input) (type check done above)
  }
  else
  {
    if (!require_dimensions(array,1) ||
        !require_native(array) || !require_contiguous(array)) SWIG_fail;
    
    SP::SiconosVector tmp;
    tmp.reset(new SiconosVector(array_size(array,0)));
    // copy : with SiconosVector based on resizable std::vector there is
    // no other way
    memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
    $1 = tmp;
  }
 }

// numpy array to SP::SimpleMatrix (here a SiconosMatrix is always a
// SimpleMatrix)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosMatrix>)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in) boost::shared_ptr<SiconosMatrix> (PyArrayObject* array=NULL, int is_new_object) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,2) ||
      !require_native(array) || !require_contiguous(array)) SWIG_fail;

  SP::SimpleMatrix tmp;
  tmp.reset(new SimpleMatrix(array_size(array,0), array_size(array,1)));
  // copy : with SimpleMatrix based on resizable std::vector
  memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
  $1 = tmp;
 }


// from C++ to python 
%template() boost::shared_ptr<SiconosVector>;
%template() boost::shared_ptr<SiconosVector>;
%template() boost::shared_ptr<BlockVector>;
%template() boost::shared_ptr<SiconosMatrix>;
%template() boost::shared_ptr<SimpleMatrix>;

%typemap(out) boost::shared_ptr<SiconosVector>
{
  npy_intp this_vector_dim[1];
  this_vector_dim[0]=$1->size();
  $result = PyArray_SimpleNewFromData(1,this_vector_dim,NPY_DOUBLE,$1->getArray());
}

%typemap(out) boost::shared_ptr<SiconosMatrix>
{
  npy_intp this_matrix_dim[2];
  this_matrix_dim[0]=$1->size(0);
  this_matrix_dim[1]=$1->size(1);
  $result = PyArray_SimpleNewFromData(2,this_matrix_dim,NPY_DOUBLE,$1->getArray());
  PyArray_UpdateFlags((PyArrayObject *)$result, NPY_FORTRAN);
}



// SiconosMatrix in api mean SimpleMatrix here
%apply (boost::shared_ptr<SiconosVector>) { (SP::SiconosVector) };
%apply (boost::shared_ptr<SiconosVector>) { (boost::shared_ptr<SiconosVector>) };
%apply (boost::shared_ptr<SiconosVector>) { (SP::SiconosVector) };

%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SimpleMatrix) };
%apply (boost::shared_ptr<SiconosMatrix>) { (boost::shared_ptr<SimpleMatrix>) };
%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SiconosMatrix) }; 
///%include "std_string.i"

%include "MBTB_PYTHON_API.hpp"
 //%include "MBTB_Body.hpp"
