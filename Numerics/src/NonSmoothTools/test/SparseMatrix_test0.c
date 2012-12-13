#include "SparseMatrix.h"


/* create an empty triplet matrix, insert 2 elements, print and free */
int main()
{
  SparseMatrix *m = cs_spalloc(0,0,0,0,1); /* coo format */
  
  int info1 = 1-cs_entry(m, 3, 4, 1.0);
  int info2 = 1-cs_entry(m, 1, 2, 2.0);

  int info3 = 1-cs_print(m, 0);

  m=cs_spfree(m);

  int info4 = 1-(m==NULL);

  return info1+info2+info3+info4;
}
