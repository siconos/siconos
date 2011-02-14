// -*- C++ -*-
// Siconos-Front-End version 3.2.0, Copyright INRIA 2005-2010.
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.	
// Siconos is a free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// Siconos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Siconos; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr 
//	
// SWIG interface for Siconos Kernel types

// numpy array to SP::SimpleVector (here a SiconosVector is always a
// dense SimpleVector)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosVector>)
{
  int res = SWIG_ConvertPtr($input, 0, SWIGTYPE_p_boost__shared_ptrT_SiconosVector_t, 0);
  int state = SWIG_CheckState(res);
  $1 = is_array($input) || PySequence_Check($input) || state;
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
    
    SP::SimpleVector tmp;
    tmp.reset(new SimpleVector(array_size(array,0)));
    // copy : with SimpleVector based on resizable std::vector there is
    // no other way
    memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
    $1 = tmp;
  }
 }


%typemap(directorout, fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> ()
{
  // directorout from KernelTypes.i
  void * swig_argp;
  int swig_res = SWIG_ConvertPtr(result,&swig_argp,SWIGTYPE_p_boost__shared_ptrT_SiconosVector_t,  0  | 0);


// strange HAlpha not initialized (LinearOSNS)
// so this is commented
//  if (!SWIG_IsOK(swig_res)) 
//  {
//    Swig::DirectorTypeMismatchException::raise(SWIG_ErrorType(SWIG_ArgError(swig_res)), "in output value of type '""SP::SiconosVector""'");
//  }

  if ((!swig_argp) || (!SWIG_IsOK(swig_res)))
  {
    return (SP::SiconosVector) c_result;
  }
  else
  {
    c_result = *(reinterpret_cast< SP::SiconosVector * >(swig_argp));
    if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SP::SiconosVector * >(swig_argp);
    return (SP::SiconosVector) c_result;
  }
}

%typemap(directorout, fragment="NumPy_Fragments") boost::shared_ptr<SiconosMatrix> ()
{
  // directorout from KernelTypes.i
  void * swig_argp;
  int swig_res = SWIG_ConvertPtr(result,&swig_argp,SWIGTYPE_p_boost__shared_ptrT_SiconosMatrix_t,  0  | 0);

// idem SiconosVector
//  if (!SWIG_IsOK(swig_res)) 
//  {
//    Swig::DirectorTypeMismatchException::raise(SWIG_ErrorType(SWIG_ArgError(swig_res)), "in output value of type '""SP::SiconosMatrix""'");
//  }

  if ((!swig_argp) || (!SWIG_IsOK(swig_res)))
  {  
    return (SP::SiconosMatrix) c_result;
  }
  else
  {
    c_result = *(reinterpret_cast< SP::SiconosMatrix * >(swig_argp));
    if (SWIG_IsNewObj(swig_res)) delete reinterpret_cast< SP::SiconosMatrix * >(swig_argp);
    return (SP::SiconosMatrix) c_result;
  }
}


// numpy array to SP::SimpleMatrix (here a SiconosMatrix is always a
// dense SimpleMatrix)
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosMatrix>)
{
  int res = SWIG_ConvertPtr($input, 0, SWIGTYPE_p_boost__shared_ptrT_SiconosMatrix_t, 0);
  int state = SWIG_CheckState(res);
  $1 = is_array($input) || PySequence_Check($input) || state;
}
%typemap(in) boost::shared_ptr<SiconosMatrix> (PyArrayObject* array=NULL, int is_new_object) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,2) ||
      !require_native(array) || !require_contiguous(array)) SWIG_fail;

  SP::SimpleMatrix tmp;
  tmp.reset(new SimpleMatrix(array_size(array,0), array_size(array,1)));
  // copy this is due to SimpleMatrix based on resizable std::vector
  memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
  $1 = tmp;
 }


// from C++ to python 
%template() boost::shared_ptr<SiconosVector>;
%template() boost::shared_ptr<SimpleVector>;
%template() boost::shared_ptr<BlockVector>;
%template() boost::shared_ptr<SiconosMatrix>;
%template() boost::shared_ptr<SimpleMatrix>;

%typemap(out) boost::shared_ptr<SiconosVector>
{
  if ($1)
  {
    npy_intp this_vector_dim[1];
    this_vector_dim[0]=$1->size();
    // warning shared_ptr counter lost here du to getArray()
    $result = PyArray_SimpleNewFromData(1,this_vector_dim,NPY_DOUBLE,$1->getArray());
  }
  else
  {
    $result = Py_None;
  }
}

%typemap(out) boost::shared_ptr<SiconosMatrix>
{
  if ($1)
  {
    npy_intp this_matrix_dim[2];
    this_matrix_dim[0]=$1->size(0);
    this_matrix_dim[1]=$1->size(1);
    // warning shared_ptr counter lost here du to getArray()
    $result = PyArray_SimpleNewFromData(2,this_matrix_dim,NPY_DOUBLE,$1->getArray());
    PyArray_UpdateFlags((PyArrayObject *)$result, NPY_FORTRAN);
  }
  else
  {
    $result = Py_None;
  }
}

// Siconos{Vector,Matrix} in api mean Simple{Vector,Matrix} here
%apply (boost::shared_ptr<SiconosVector>) { (SP::SimpleVector) };
%apply (boost::shared_ptr<SiconosVector>) { (boost::shared_ptr<SimpleVector>) };
%apply (boost::shared_ptr<SiconosVector>) { (SP::SiconosVector) };

%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SimpleMatrix) };
%apply (boost::shared_ptr<SiconosMatrix>) { (boost::shared_ptr<SimpleMatrix>) };
%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SiconosMatrix) }; 

