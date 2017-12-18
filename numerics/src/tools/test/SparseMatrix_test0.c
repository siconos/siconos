
#include "SparseMatrix_internal.h"

/* create an empty triplet matrix, insert 2 elements, print and free */
int main()
{
  CSparseMatrix *m = cs_spalloc(0,0,0,0,1); /* coo format */

  CS_INT info1 = 1-cs_entry(m, 3, 4, 1.0);
  CS_INT info2 = 1-cs_entry(m, 1, 2, 2.0);

  CS_INT info3 = 1-cs_print(m, 0);

  m=cs_spfree(m);

  CS_INT info4 = 1-(m==NULL);

  return (int)(info1+info2+info3+info4);
}
