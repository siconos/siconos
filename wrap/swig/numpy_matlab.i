

%typemap(in) (double IN_ARRAY1[ANY]) (SN_ARRAY_TYPE* array=NULL, int is_new_object=0) {
  array = obj_to_sn_vector($input, &is_new_object);
  if(array && CHECK_ARRAY_VECTOR(array)) SWIG_fail;
  $1 = ($1_ltype) array_data(array);
}

/* Typemap suite for (DATA_TYPE ARGOUT_ARRAY1[ANY])
 */
%typemap(in,numinputs=0)
  (double ARGOUT_ARRAY1[ANY])
  (SN_ARRAY_TYPE* array = NULL)
{
  C_to_target_lang1_alloc($1, array, $1_dim0, SWIG_fail)
}
%typemap(argout)
  (double ARGOUT_ARRAY1[ANY])
{
  $result = SWIG_AppendOutput($result,(SN_OBJ_TYPE*)array$argnum);
}

/* Typemap suite for (DIM_TYPE DIM1, DATA_TYPE* INPLACE_ARRAY1)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY)
  (int DIM1, double* INPLACE_ARRAY1)
{
  $1 = mxIsDouble($input)
}
%typemap(in) (int DIM1, double* INPLACE_ARRAY1) (SN_ARRAY_TYPE* array=NULL, int i=0)
{ 
  array = obj_to_sn_vector($input, &i);
  if (CHECK_ARRAY_VECTOR(array)) SWIG_fail;
  $1 = 1;
  for (i=0; i < array_numdims(array); ++i) $1 *= array_size(array,i);
  $2 = (double*) array_data(array);
}

/* Typemap suite for (int DIM1, int DIM2, double* IN_FARRAY2)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (int DIM1, int DIM2, double* IN_FARRAY2)
{
  $1 = is_array($input) || PySequence_Check($input);
}
%typemap(in,
         fragment="NumPy_Fragments")
  (int DIM1, int DIM2, double* IN_FARRAY2)
  (SN_ARRAY_TYPE* array=NULL, int is_new_object=0)
{

  array = obj_to_sn_vector($input, &is_new_object);

  if (CHECK_ARRAY_MATRIX(array)) SWIG_fail;

  $1 = (int) array_size(array,0);
  $2 = (int) array_size(array,1);
  $3 = (double*) array_data(array);
}

