// vectors of size problem_size from given *Problem as first input
// no conversion => inout array XXX FIX issue here
%typemap(in) (double *z) (PyArrayObject* array=NULL, int is_new_object=0) {
  //typemap(in) (double *z) (PyArrayObject* array=NULL, int is_new_object=0)

  array = obj_to_array_allow_conversion($input, NPY_DOUBLE, &is_new_object);

  if (!array || !require_native(array) || !require_contiguous(array) || !require_fortran(array)) SWIG_fail;

  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *z)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// list of matrices problemSizex3
%typemap(in) (double *blocklist3x3) (PyArrayObject* array=NULL, int is_new_object=0) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  npy_intp array_len[2] = {0,0};

  if (! *p_problem_size1)
  {
    if (array_numdims(array) > 1)
    {
      *p_problem_size1 = fmax(array_size(array,0), array_size(array,1)) / 3;
    }
    else if (array_numdims(array) == 1)
    {
      *p_problem_size1 = array_size(array,0) / 3;
    }
    else
      SWIG_fail;

    if (*p_problem_size1 % 3 != 0) SWIG_fail;

    if (*p_problem_size1 / 3 == 0) SWIG_fail;

    number_of_contacts1 = *p_problem_size1 / 3;

  }

  assert (*p_problem_size1);


  if (array_numdims(array) > 1)
  {
    array_len[0] = *p_problem_size1 * 3;
    array_len[1] = 1;
  }
  else
  {
    array_len[0] = *p_problem_size1 * 3;
  }

  if (!array
      || !require_native(array) || !require_contiguous(array)
      || !require_size(array, array_len, array_numdims(array))) SWIG_fail;

  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blocklist3x3)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}


// matrices problemSizexproblemSize
%typemap(in) (double *blockarray3x3) (PyArrayObject* array=NULL, int is_new_object=0) {

  array = obj_to_array_fortran_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  npy_intp array_len[2] = {0,0};

  if (! *p_problem_size1)
  {
    if (array_numdims(array) > 1)
    {
      *p_problem_size1 = array_size(array,0); // assume square matrix
    }
    else if (array_numdims(array) == 1)
    {
      *p_problem_size1 = isqrt(array_size(array,0));
    }
    else
      SWIG_fail;


    if (*p_problem_size1 % 3 != 0) SWIG_fail;

    if (*p_problem_size1 / 3 == 0) SWIG_fail;

    number_of_contacts1 = *p_problem_size1 / 3;

  }

  assert (*p_problem_size1);

  if (array_numdims(array) > 1)
  {
    array_len[0] = *p_problem_size1;
    array_len[1] = *p_problem_size1;
  }
  else
  {
    array_len[0] = *p_problem_size1 * *p_problem_size1;
  }

  if (!array
      || !require_native(array) || !require_fortran(array)
      || !require_size(array, array_len, array_numdims(array))) SWIG_fail;

  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blockarray3x3)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// vectors of size problem_size from problemSize as first input
%typemap(in) (double *blocklist3) (PyArrayObject* array=NULL, int is_new_object=0) {

  array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE,&is_new_object);

  npy_intp array_len[2] = {0,0};

  if (! *p_problem_size1)
  {
    if (array_numdims(array) == 1)
    {
      *p_problem_size1 = array_size(array,0);
    }
    else if (array_numdims(array) > 1)
    {
      *p_problem_size1 = fmax(array_size(array,0), array_size(array,1));
    }

    if (*p_problem_size1 % 3 != 0) SWIG_fail;

    if (*p_problem_size1 / 3 == 0) SWIG_fail;


    number_of_contacts1 = *p_problem_size1 / 3;

  }

  assert (*p_problem_size1);

  if (array_numdims(array) == 1)
  {
    array_len[0] = *p_problem_size1;
  }
  else
  {
    array_len[0] = *p_problem_size1;
    array_len[1] = 1;
  }

  if (!array
      || !require_native(array) || !require_contiguous(array)
      || !require_size(array, array_len, array_numdims(array))) SWIG_fail;

  $1 = (double *) array_data(array);

 }

%typemap(freearg) (double *blocklist3)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

// vectors of size problem_size

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *output_blocklist3) (PyObject* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

// 2 : check must be done after in
%typemap(check) (double *output_blocklist3)
{
  if (*p_problem_size1)
  {

    npy_intp dims[2] = { *p_problem_size1, 1};

    array$argnum = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (!array$argnum) SWIG_fail;
    $1 = ($1_ltype) array_data(array$argnum);
  }

}

// 3 : return arg
%typemap(argout) (double *output_blocklist3)
{
  if (*p_problem_size1)
  {
     $result = SWIG_Python_AppendOutput($result,(PyObject *)array$argnum);
  }

}

// 3x3 matrices

// 1 : numinputs=0 mandatory to avoid arg
%typemap(in, numinputs=0) (double *output_blocklist3x3) (PyObject* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

%typemap(in, numinputs=0) (double *output_blockarray3x3) (PyObject* array=NULL)
{
    // %typemap(in, numinputs=0)
    // we cannot get problem_size here as numinputs=0 => before
    // numinputs=1, how can we change this ??
}

// 2 : check must be done after in
%typemap(check) (double *output_blocklist3x3)
{
  if (*p_problem_size1)
  {

    npy_intp dims[2] = { *p_problem_size1 * 3, 1 };

    array$argnum = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    // block list : require_fortran useless?
    if (!array$argnum) SWIG_fail;
    PyArrayObject *array = (PyArrayObject*) array$argnum;
    if (!array || !require_fortran(array)) SWIG_fail;
    $1 = ($1_ltype) array_data(array);
  }

}

%typemap(check) (double *output_blockarray3x3)
{
  if (*p_problem_size1)
  {

    npy_intp dims[2] = { *p_problem_size1, *p_problem_size1};

    array$argnum = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    PyArrayObject *array = (PyArrayObject*) array$argnum;
    if (!array || !require_fortran(array)) SWIG_fail;
    $1 = ($1_ltype) array_data(array);
  }

}

// 3 : return arg
%typemap(argout) (double *output_blocklist3x3)
{
  if (*p_problem_size1)
  {
    $result = SWIG_Python_AppendOutput($result,(PyObject *)array$argnum);
  }

}

%typemap(argout) (double *output_blockarray3x3)
{
  if (*p_problem_size1)
  {
    $result = SWIG_Python_AppendOutput($result,(PyObject *)array$argnum);
  }
}

// vectors of size numberOfContacts
%typemap(memberin) (double *mu) {
  // Still some dark magic :( --xhub
  if (arg1->numberOfContacts <= 0)
  {
    PyErr_SetString(PyExc_RuntimeError, "numberOfContacts is not set, it sould be done first!");
    SWIG_fail;
  }

  if (arg1->numberOfContacts !=  array_size(array2, 0))
  {
    char msg[1024];
    snprintf(msg, sizeof(msg), "Size of mu is %ld, but the number of contacts is %d! Both should be equal!\n", array_size(array2, 0), arg1->numberOfContacts);
    PyErr_SetString(PyExc_RuntimeError, msg);
    SWIG_fail;
  }

  if (!$1) { $1 = (double*)malloc(arg1->numberOfContacts * sizeof(double)); }
  memcpy($1, $input, arg1->numberOfContacts * sizeof(double));

 }

// vectors of size M
%typemap(memberin) (double *q) {
  // Still some dark magic :( --xhub
 char msg[1024];
  assert(arg1);
  if (!arg1->M)
  {
    PyErr_SetString(PyExc_RuntimeError, "M is not initialized, it sould be done first!");
    SWIG_fail;
  }

  int size = arg1->M->size0;
  if (size !=  array_size(array2, 0))
  {
    snprintf(msg, sizeof(msg), "Size of q is %ld, but the size of M is %d! Both should be equal!\n", array_size(array2, 0), size);
    PyErr_SetString(PyExc_RuntimeError, msg);
    SWIG_fail;
  }

  if (!$1) { $1 = (double*)malloc(size * sizeof(double)); }
  memcpy($1, $input, size * sizeof(double));

 }

%typemap(in, numinputs=0) (double *output3) 
{
}

%typemap(argout) (double *output3)
{
  $result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble($1[0]));
  $result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble($1[1]));
  $result = SWIG_Python_AppendOutput($result, PyFloat_FromDouble($1[2]));
}

// other names that must be transformed this way
%apply (double *z) { (double *r) };

%apply (double *z) { (double *u) };

%apply (double *z) { (double *w) };

%apply (double *z) { (double *x) };

%apply (double *z) { (double *F) };

%apply (double *z) { (double *Fmcp) };

%apply (double *z) { (double *zlem) };

%apply (double *z) { (double *wlem) };

%apply (double *z) { (double *reaction) };

%apply (double *z) { (double *velocity) };

%apply (double *z) { (double *globalVelocity) };

%apply (double *z) { (double *mu) };

%apply (double *z) { (double *q) };

%apply (double *z) { (double *b) };

%apply (double *blocklist3) { (double *vect3D) };

%apply (double *blocklist3) { (double *reaction3D) };

%apply (double *blocklist3) { (double *velocity3D) };

%apply (double *blocklist3) { (double *rho3D) };

%apply (double *blocklist3) { (double *blocklist3_1) };

%apply (double *blocklist3) { (double *blocklist3_2) };

%apply (double *blocklist3) { (double *blocklist3_3) };

%apply (double *blocklist3x3) { (double *blocklist3x3_1) };

%apply (double *blocklist3x3) { (double *blocklist3x3_2) };

%apply (double *blocklist3x3) { (double *blocklist3x3_3) };

%apply (double *output_blocklist3) { (double *output_blocklist3_1) };

%apply (double *output_blocklist3) { (double *output_blocklist3_2) };

%apply (double *output_blocklist3) { (double *output_blocklist3_3) };

%apply (double *output_blocklist3) { (double *output_blocklist3_4) };

%apply (double *output_blocklist3x3) { (double *output_blocklist3x3_1) };

%apply (double *output_blocklist3x3) { (double *output_blocklist3x3_2) };
