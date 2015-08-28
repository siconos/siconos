
 %typemap(in) (LinearComplementarityProblem*) (npy_intp problem_size) {
   void *lcp;
   int res = SWIG_ConvertPtr($input, &lcp,SWIGTYPE_p_LinearComplementarityProblem, 0 |  0 );
   if (!SWIG_IsOK(res)) SWIG_fail;

   problem_size=((LinearComplementarityProblem *) lcp)->size;

   $1 = (LinearComplementarityProblem *) lcp;
 }

 %typemap(in) (MixedLinearComplementarityProblem*) (npy_intp mlcproblem_size) {
   void *mlcp;
   int res = SWIG_ConvertPtr($input, &mlcp,SWIGTYPE_p_MixedLinearComplementarityProblem, 0 |  0 );
   if (!SWIG_IsOK(res)) SWIG_fail;

   mlcproblem_size=((MixedLinearComplementarityProblem *) mlcp)->n +((MixedLinearComplementarityProblem *) mlcp)->m;

   $1 = (MixedLinearComplementarityProblem *) mlcp;
 }

 // problemSize given as first arg => set by first numpy array length
 // in remaining args
 %typemap(in, numinputs=0) 
   (unsigned int problemSize) 
   (unsigned int *p_problem_size, npy_intp number_of_contacts)
 {
   // the first array length sets problemSize
   p_problem_size = &$1;
   *p_problem_size = 0;
   number_of_contacts = 0;  
 }

%typemap(in) 
  (FrictionContactProblem*) 
  (npy_intp problem_size, npy_intp problem_dimension, npy_intp number_of_contacts) 
{
  void *fcp;
  int res = SWIG_ConvertPtr($input, &fcp,SWIGTYPE_p_FrictionContactProblem, 0 |  0 );
  if (!SWIG_IsOK(res)) SWIG_fail;

  problem_dimension=((FrictionContactProblem *) fcp)->dimension;
  number_of_contacts=((FrictionContactProblem *) fcp)->numberOfContacts;
  problem_size=((FrictionContactProblem *) fcp)->numberOfContacts * problem_dimension;

  $1 = (FrictionContactProblem*) fcp;
}

%typemap(in) 
  (GlobalFrictionContactProblem*) 
  (npy_intp problem_size, npy_intp problem_dimension, npy_intp number_of_contacts) 
{
  void *fcp;
  int res = SWIG_ConvertPtr($input, &fcp,SWIGTYPE_p_GlobalFrictionContactProblem, 0 |  0 );
  if (!SWIG_IsOK(res)) SWIG_fail;

  problem_dimension=((GlobalFrictionContactProblem *) fcp)->dimension;
  number_of_contacts=((GlobalFrictionContactProblem *) fcp)->numberOfContacts;
  problem_size=((GlobalFrictionContactProblem *) fcp)->numberOfContacts * problem_dimension;

  $1 = (GlobalFrictionContactProblem*) fcp;
}


