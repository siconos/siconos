// -*- C++ -*-
// Siconos-Front-End , Copyright INRIA 2005-2011.
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



// check on input : a numpy array or a SiconosVector
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosVector>)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (boost::shared_ptr<SiconosVector>)
  int res = SWIG_ConvertPtr($input, 0, SWIGTYPE_p_boost__shared_ptrT_SiconosVector_t, 0);
  int state = SWIG_CheckState(res);
  $1 = is_array($input) || PySequence_Check($input) || state;
}

// check on input : a numpy array or a SimpleVector
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SimpleVector>)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (boost::shared_ptr<SimpleVector>)
  int res = SWIG_ConvertPtr($input, 0, SWIGTYPE_p_boost__shared_ptrT_SimpleVector_t, 0);
  int state = SWIG_CheckState(res);
  $1 = is_array($input) || PySequence_Check($input) || state;
}

// check on input : a numpy array or a SiconosMatrix
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SiconosMatrix>)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (boost::shared_ptr<SiconosMatrix>)
  int res = SWIG_ConvertPtr($input, 0, SWIGTYPE_p_boost__shared_ptrT_SiconosMatrix_t, 0);
  int state = SWIG_CheckState(res);
  $1 = is_array($input) || PySequence_Check($input) || state;
}

// check on input : a numpy array or a SimpleMatrix
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
(boost::shared_ptr<SimpleMatrix>)
{
  // %typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (boost::shared_ptr<SimpleMatrix>)
  int res = SWIG_ConvertPtr($input, 0, SWIGTYPE_p_boost__shared_ptrT_SimpleMatrix_t, 0);
  int state = SWIG_CheckState(res);
  $1 = is_array($input) || PySequence_Check($input) || state;
}



// numpy or SP::SiconosVector on input -> SP::SiconosVector
%typemap(in,fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> (PyArrayObject* array=NULL, int is_new_object)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %boost::shared_ptr<SiconosVector> (PyArrayObject* array=NULL, int
  // %is_new_object)
  void *argp1=0;
  int res1=0;
  int newmem = 0;
  boost::shared_ptr< SiconosVector > tempshared1 ;
  boost::shared_ptr< SiconosVector > *smartarg1 = 0 ;
 
  // try a conversion from a SiconosVector
  res1 = SWIG_ConvertPtrAndOwn($input, &argp1, SWIGTYPE_p_boost__shared_ptrT_SiconosVector_t, 0 |  0 , &newmem);
  if (SWIG_IsOK(res1)) 
  {
    if (newmem & SWIG_CAST_NEW_MEMORY) 
    {
      // taken from generated code, we assume we should never get here
      assert(false);
      tempshared1 = *reinterpret_cast< boost::shared_ptr< SiconosVector > * >(argp1);
      delete reinterpret_cast< boost::shared_ptr< SiconosVector > * >(argp1);
    } 
    else {
      smartarg1 = reinterpret_cast< boost::shared_ptr< SiconosVector > * >(argp1);
      $1 = *smartarg1;
    }
  }
  else
  {
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

    if (!array)
    {
      void *argp;
      SWIG_fail; // not implemented : $1 = type_conv($input) (type check done above)
    }
    else
    {
      if (!require_dimensions(array,1) ||
          !require_native(array) || !require_fortran(array)) SWIG_fail;
      
      SP::SimpleVector tmp;
      tmp.reset(new SimpleVector(array_size(array,0)));
      // copy : with SimpleVector based on resizable std::vector there is
      // no other way
      memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
      $1 = tmp;
    }
  }
}

// numpy or SP::SimpleVector on input -> SP::SimpleVector
%typemap(in,fragment="NumPy_Fragments") boost::shared_ptr<SimpleVector> (PyArrayObject* array=NULL, int is_new_object)
{
  // %typemap(in,fragment="NumPy_Fragments")
  // %boost::shared_ptr<SimpleVector> (PyArrayObject* array=NULL, int
  // %is_new_object)
  void *argp1=0;
  int res1=0;
  int newmem = 0;
  boost::shared_ptr< SimpleVector > tempshared1 ;
  boost::shared_ptr< SimpleVector > *smartarg1 = 0 ;
 
  // try a conversion from a SimpleVector
  res1 = SWIG_ConvertPtrAndOwn($input, &argp1, SWIGTYPE_p_boost__shared_ptrT_SimpleVector_t, 0 |  0 , &newmem);
  if (SWIG_IsOK(res1)) 
  {
    if (newmem & SWIG_CAST_NEW_MEMORY) 
    {
      // taken from generated code, we assume we should never get here
      assert(false);
      tempshared1 = *reinterpret_cast< boost::shared_ptr< SimpleVector > * >(argp1);
      delete reinterpret_cast< boost::shared_ptr< SimpleVector > * >(argp1);
    } 
    else {
      smartarg1 = reinterpret_cast< boost::shared_ptr< SimpleVector > * >(argp1);
      $1 = *smartarg1;
    }
  }
  else
  {
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

    if (!array)
    {
      void *argp;
      SWIG_fail; // not implemented : $1 = type_conv($input) (type check done above)
    }
    else
    {
      if (!require_dimensions(array,1) ||
          !require_native(array) || !require_fortran(array)) SWIG_fail;
      
      SP::SimpleVector tmp;
      tmp.reset(new SimpleVector(array_size(array,0)));
      // copy : with SimpleVector based on resizable std::vector there is
      // no other way
      memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
      $1 = tmp;
    }
  }
}

// director input : SP::SiconosVector -> numpy
%typemap(directorin, fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> ()
{
  // %typemap(directorin, fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> ()
  if($1_name)
  {
    if (ask<IsDense>(*$1_name))
    {
      npy_intp this_vector_dim[1];
      this_vector_dim[0]=$1_name->size();
      // warning shared_ptr counter lost here du to getArray()
      $input = PyArray_SimpleNewFromData(1,this_vector_dim,NPY_DOUBLE,$1_name->getArray());
    }
    else
    {
      // not a dense vector : no conversion
      $input = SWIG_NewPointerObj(SWIG_as_voidptr(&$1_name), SWIGTYPE_p_boost__shared_ptrT_SiconosVector_t,  0 );  
    }
  }
  else
  {
    Py_INCREF(Py_None);
    $input = Py_None;
  }
}

// director input : SP::SimpleVector -> numpy
%typemap(directorin, fragment="NumPy_Fragments") boost::shared_ptr<SimpleVector> ()
{
  // %typemap(directorin, fragment="NumPy_Fragments") boost::shared_ptr<SimpleVector> ()
  if($1_name)
  {
    if (ask<IsDense>(*$1_name))
    {
      npy_intp this_vector_dim[1];
      this_vector_dim[0]=$1_name->size();
      // warning shared_ptr counter lost here du to getArray()
      $input = PyArray_SimpleNewFromData(1,this_vector_dim,NPY_DOUBLE,$1_name->getArray());
    }
    else
    {
      // not a dense vector : no conversion
      $input = SWIG_NewPointerObj(SWIG_as_voidptr(&$1_name), SWIGTYPE_p_boost__shared_ptrT_SimpleVector_t,  0 );  
    }
  }
  else
  {
    Py_INCREF(Py_None);
    $input = Py_None;
  }
}

// director input : SP::SiconosMatrix -> numpy
%typemap(directorin, fragment="NumPy_Fragments") boost::shared_ptr<SiconosMatrix> ()
{
  // %typemap(directorin, fragment="NumPy_Fragments") boost::shared_ptr<SiconosMatrix> ()
  if ($1_name)
  {
    if (($1_name).getNum() == 1)
    {
      npy_intp this_matrix_dim[2];
      this_matrix_dim[0]=$1->size(0);
      this_matrix_dim[1]=$1->size(1);
      // warning shared_ptr counter lost here du to getArray()
      $input = PyArray_SimpleNewFromData(2,this_matrix_dim,NPY_DOUBLE,$1->getArray());
      PyArray_UpdateFlags((PyArrayObject *)$input, NPY_FORTRAN);
    }
    else
    {
      // not a dense matrix : no conversion
      $input = SWIG_NewPointerObj(SWIG_as_voidptr(&$1_name), SWIGTYPE_p_boost__shared_ptrT_SiconosMatrix_t,  0 );  
    }

  }
  else
  {
    Py_INCREF(Py_None);
    $input = Py_None;
  }
}

// director input : SP::SimpleMatrix -> numpy
%typemap(directorin, fragment="NumPy_Fragments") boost::shared_ptr<SimpleMatrix> ()
{
  // %typemap(directorin, fragment="NumPy_Fragments") boost::shared_ptr<SimpleMatrix> ()
  if ($1_name)
  {
    if (($1_name).getNum() == 1)
    {
      npy_intp this_matrix_dim[2];
      this_matrix_dim[0]=$1->size(0);
      this_matrix_dim[1]=$1->size(1);
      // warning shared_ptr counter lost here du to getArray()
      $input = PyArray_SimpleNewFromData(2,this_matrix_dim,NPY_DOUBLE,$1->getArray());
      PyArray_UpdateFlags((PyArrayObject *)$input, NPY_FORTRAN);
    }
    else
    {
      // not a dense matrix : no conversion
      $input = SWIG_NewPointerObj(SWIG_as_voidptr(&$1_name), SWIGTYPE_p_boost__shared_ptrT_SimpleMatrix_t,  0 );  
    }
  }
  else
  {
    Py_INCREF(Py_None);
    $input = Py_None;
  }
}

// director output : PyObject -> SP::SiconosVector 
%typemap(directorout, fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> ()
{
  // %typemap(directorout, fragment="NumPy_Fragments") boost::shared_ptr<SiconosVector> ()
  void * swig_argp;

  // try a conversion from SP::SiconosVector
  int swig_res = SWIG_ConvertPtr(result,&swig_argp,SWIGTYPE_p_boost__shared_ptrT_SiconosVector_t,  0  | 0);

  if (!SWIG_IsOK(swig_res))
  {
    // try a conversion from numpy
    PyArrayObject* array = NULL;
    int is_new_object;
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);
    if (!require_dimensions(array,1) ||
        !require_native(array) || !require_fortran(array)) throw Swig::DirectorMethodException();
    
    SP::SimpleVector tmp;
    tmp.reset(new SimpleVector(array_size(array,0)));
    // copy : with SimpleVector based on resizable std::vector there is
    // no other way
    memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*sizeof(double));
    return tmp;
  }

  if (!swig_argp)
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

// director output : PyObject -> SP::SimpleVector 
%typemap(directorout, fragment="NumPy_Fragments") boost::shared_ptr<SiconosMatrix> ()
{
  // %typemap(directorout, fragment="NumPy_Fragments") boost::shared_ptr<SiconosMatrix> ()
  void * swig_argp;
  int swig_res = SWIG_ConvertPtr(result,&swig_argp,SWIGTYPE_p_boost__shared_ptrT_SiconosMatrix_t,  0  | 0);

  if (!SWIG_IsOK(swig_res))
  {
    // try a conversion from numpy
    PyArrayObject* array = NULL;
    int is_new_object;
    array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);
    if (!require_dimensions(array,2) ||
        !require_native(array) || !require_fortran(array)) throw Swig::DirectorMethodException();
    

    SP::SimpleMatrix tmp;
    tmp.reset(new SimpleMatrix(array_size(array,0), array_size(array,1)));
    // copy this is due to SimpleMatrix based on resizable std::vector
    memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
    return tmp;
  }

  if (!swig_argp)
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


%typemap(in) boost::shared_ptr<SiconosMatrix> (PyArrayObject* array=NULL, int is_new_object) {

  // %typemap(in) boost::shared_ptr<SiconosMatrix> (PyArrayObject* array=NULL, int is_new_object)
  array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,2) ||
      !require_native(array) || !require_fortran(array)) SWIG_fail;

  SP::SimpleMatrix tmp;
  tmp.reset(new SimpleMatrix(array_size(array,0), array_size(array,1)));
  // copy this is due to SimpleMatrix based on resizable std::vector
  memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
  $1 = tmp;
 }

%typemap(in) boost::shared_ptr<SimpleMatrix> (PyArrayObject* array=NULL, int is_new_object) {

  // %typemap(in) boost::shared_ptr<SiconosMatrix> (PyArrayObject* array=NULL, int is_new_object)
  array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  if (!array || !require_dimensions(array,2) ||
      !require_native(array) || !require_fortran(array)) SWIG_fail;

  SP::SimpleMatrix tmp;
  tmp.reset(new SimpleMatrix(array_size(array,0), array_size(array,1)));
  // copy this is due to SimpleMatrix based on resizable std::vector
  memcpy(&*tmp->getArray(),array_data(array),array_size(array,0)*array_size(array,1)*sizeof(double));
  $1 = tmp;
 }

%inline %{
  // SWIG_DIRECTOR_CAST cannot be use on non polymorphic data type
  // (i.e. classes without at least one virtual method) so we have to
  // switch at compile time. This is done with boost::mpl::eval_if

  template<typename T>
    struct DirectorCast
  {
    typedef DirectorCast type;
    Swig::Director* value(T& p)
    {
      return SWIG_DIRECTOR_CAST(p);
    }
  };
  
  template<typename T>
    struct DirectorNoCast
  {
    typedef DirectorNoCast type;
    Swig::Director* value(T& p)
    {
      return 0;
    }
  };
%}  
  


%typemap(out) boost::shared_ptr<SiconosVector>
{
  // %typemap(out) boost::shared_ptr<SiconosVector>

  // compile time test to reproduce swig director test. swig does not
  // seem to provides facilities to customize this
  // arg1 is the class instantiation (swig 2.0) => no way to get this as a swig
  // variable.
  typedef BOOST_TYPEOF(arg1) self_type;

  typedef boost::mpl::eval_if<boost::is_polymorphic<self_type >,
    DirectorCast<self_type >,
    DirectorNoCast<self_type > >::type CastMaybe;
  
  CastMaybe cast_maybe;
  
  Swig::Director* l_director = cast_maybe.value(arg1);
  bool l_upcall = (l_director && (l_director->swig_get_self()==obj0));

  // call from director?
  if (l_upcall)
  {
    // result from C++ method, return the pointer
    $result = SWIG_NewPointerObj(SWIG_as_voidptr(&$1), SWIGTYPE_p_boost__shared_ptrT_SiconosVector_t,  0 );
  }
  // call from python : return numpy from SiconosVector
  else
  {
    if ($1)
    {
      // /!\ need check for a dense vector!

      npy_intp this_vector_dim[1];
      this_vector_dim[0]=$1->size();
      // warning shared_ptr counter lost here du to getArray()
      $result = PyArray_SimpleNewFromData(1,this_vector_dim,NPY_DOUBLE,$1->getArray());
    }
    else
    {
      Py_INCREF(Py_None);
      $result = Py_None;
    }
  }
}

%typemap(out) boost::shared_ptr<SimpleVector> (bool upcall=false)
{
  // %typemap(out) boost::shared_ptr<SimpleVector>

  // compile time test to reproduce swig director test. swig does not
  // seem to provides facilities to customize this
  // arg1 is the class instantiation (swig 2.0) => no way to get this as a swig
  // variable.
  typedef BOOST_TYPEOF(arg1) self_type;

  typedef boost::mpl::eval_if<boost::is_polymorphic<self_type >,
    DirectorCast<self_type >,
    DirectorNoCast<self_type > >::type CastMaybe;
  
  CastMaybe cast_maybe;
  
  Swig::Director* l_director = cast_maybe.value(arg1);
  bool l_upcall = (l_director && (l_director->swig_get_self()==obj0));

  // call from director?
  if (l_upcall)
  {
    // result from C++ method, return the pointer
    $result = SWIG_NewPointerObj(SWIG_as_voidptr(&$1), SWIGTYPE_p_boost__shared_ptrT_SimpleVector_t,  0 );
  }
  // call from python : return numpy from SimpleVector
  else
  {
    if ($1)
    {
      // /!\ need check for a dense vector!

      npy_intp this_vector_dim[1];
      this_vector_dim[0]=$1->size();
      // warning shared_ptr counter lost here du to getArray()
      $result = PyArray_SimpleNewFromData(1,this_vector_dim,NPY_DOUBLE,$1->getArray());
    }
    else
    {
      Py_INCREF(Py_None);
      $result = Py_None;
    }
  }
}


//PyArray_UpdateFlags does not seem to have any effect
//>>> r = K.FirstOrderLinearTIR()
//>>> r.setCPtr([[1,2,3],[4,5,6]])
//>>> C=r.C()
//>>> C.flags
//  C_CONTIGUOUS : True
//  F_CONTIGUOUS : False            <---- !!!
//  OWNDATA : False
//  WRITEABLE : True
//  ALIGNED : True
//  UPDATEIFCOPY : False
//
//
// with this macro : ok
#define FPyArray_SimpleNewFromData(nd, dims, typenum, data)             \
  PyArray_New(&PyArray_Type, nd, dims, typenum, NULL,                   \
              data, 0, NPY_FARRAY, NULL)


%typemap(out) boost::shared_ptr<SiconosMatrix> (bool upcall=false)
{
  // %typemap(out) boost::shared_ptr<SiconosMatrix>

  // compile time test to reproduce swig director test. swig does not
  // seem to provides facilities to customize this
  // arg1 is the class instantiation (swig 2.0) => no way to get this as a swig
  // variable.
  typedef BOOST_TYPEOF(arg1) self_type;

  typedef boost::mpl::eval_if<boost::is_polymorphic<self_type >,
    DirectorCast<self_type >,
    DirectorNoCast<self_type > >::type CastMaybe;
  
  CastMaybe cast_maybe;
  
  Swig::Director* l_director = cast_maybe.value(arg1);
  bool l_upcall = (l_director && (l_director->swig_get_self()==obj0));

  // call from director?
  if (l_upcall)
  {
    // result from C++ method, return the pointer
    $result = SWIG_NewPointerObj(SWIG_as_voidptr(&$1), SWIGTYPE_p_boost__shared_ptrT_SiconosMatrix_t,  0 );
  }
  // call from python : return numpy from SiconosMatrix
  else
  {
    if ($1)
    {
      // /!\ need check for a dense matrix!

      npy_intp this_matrix_dim[2];
      this_matrix_dim[0]=$1->size(0);
      this_matrix_dim[1]=$1->size(1);
      // warning shared_ptr counter lost here du to getArray()
      $result = FPyArray_SimpleNewFromData(2,this_matrix_dim,NPY_DOUBLE,$1->getArray());

      // see comments about PyArray_UpdateFlags above
      PyArray_UpdateFlags((PyArrayObject *)$result, (PyArray_FLAGS((PyArrayObject *)$result))|NPY_CONTIGUOUS|NPY_FORTRAN);
    }
    else
    {
      Py_INCREF(Py_None);
      $result = Py_None;
    }
  }
}

%typemap(out) boost::shared_ptr<SimpleMatrix> (bool upcall=false)
{
  // %typemap(out) boost::shared_ptr<SimpleMatrix>

  // compile time test to reproduce swig director test. swig does not
  // seem to provides facilities to customize this
  // arg1 is the class instantiation (swig 2.0) => no way to get this as a swig
  // variable.
  typedef BOOST_TYPEOF(arg1) self_type;

  typedef boost::mpl::eval_if<boost::is_polymorphic<self_type >,
    DirectorCast<self_type >,
    DirectorNoCast<self_type > >::type CastMaybe;
  
  CastMaybe cast_maybe;
  
  Swig::Director* l_director = cast_maybe.value(arg1);
  bool l_upcall = (l_director && (l_director->swig_get_self()==obj0));

  // call from director?
  if (l_upcall)
  {
    // result from C++ method, return the pointer
    $result = SWIG_NewPointerObj(SWIG_as_voidptr(&$1), SWIGTYPE_p_boost__shared_ptrT_SimpleMatrix_t,  0 );
  }
  // call from python : return numpy from SiconosMatrix
  else
  {
    if ($1)
    {
      // /!\ need check for a dense matrix!

      npy_intp this_matrix_dim[2];
      this_matrix_dim[0]=$1->size(0);
      this_matrix_dim[1]=$1->size(1);
      // warning shared_ptr counter lost here du to getArray()
      $result = FPyArray_SimpleNewFromData(2,this_matrix_dim,NPY_DOUBLE,$1->getArray());

      // see comments about PyArray_UpdateFlags above
      PyArray_UpdateFlags((PyArrayObject *)$result, (PyArray_FLAGS((PyArrayObject *)$result))|NPY_CONTIGUOUS|NPY_FORTRAN);
    }
    else
    {
      Py_INCREF(Py_None);
      $result = Py_None;
    }
  }
}

// check on input : a python sequence
%typecheck(SWIG_TYPECHECK_INTEGER)
(boost::shared_ptr<std::vector<unsigned int> >) 
{
  // %typecheck(boost::shared_ptr<std::vector<unsigned int> >, precedence=SWIG_TYPECHECK_INTEGER))
  PySequence_Check($input);
}


// int sequence => std::vector<unsigned int>
%typemap(in,fragment="NumPy_Fragments") boost::shared_ptr<std::vector<unsigned int> > (boost::shared_ptr<std::vector<unsigned int> > temp) 
{
  temp.reset(new std::vector<unsigned int>());
  if (!sequenceToUnsignedIntVector($input, temp))
  {
    SWIG_fail;
  }
  $1 = temp; // temp deallocation is done at object destruction
             // thanks to shared ptr ref counting
}
  
// needed?
// from C++ to python 
%template() boost::shared_ptr<SiconosVector>;
%template() boost::shared_ptr<SimpleVector>;
%template() boost::shared_ptr<BlockVector>;
%template() boost::shared_ptr<SiconosMatrix>;
%template() boost::shared_ptr<SimpleMatrix>;
%template() boost::shared_ptr<std::vector<unsigned int> >;


%apply (boost::shared_ptr<SiconosVector>) { (SP::SiconosVector) };
%apply (boost::shared_ptr<SimpleVector>) { (SP::SimpleVector) };

%apply (boost::shared_ptr<SiconosMatrix>) { (SP::SiconosMatrix) };
%apply (boost::shared_ptr<SimpleMatrix>) { (SP::SimpleMatrix) };

%apply (boost::shared_ptr<std::vector<unsigned int> >) { (SP::UnsignedIntVector) };

