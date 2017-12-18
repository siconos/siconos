 // problemSize given as first arg => set by first numpy array length
 // in remaining args
 %typemap(in, numinputs=0) 
   (unsigned int problemSize) 
   (unsigned int *p_problem_size, SN_ARRAY_INT_TYPE number_of_contacts)
 {
   // the first array length sets problemSize
   p_problem_size = &$1;
   *p_problem_size = 0;
   number_of_contacts = 0;  
 }
