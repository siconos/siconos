#include "op3x3.h"


/** print a matrix
 * \param double* a
 */
#include <stdio.h>
void print3x3(double* mat)
{
  SET3X3(mat);

  printf("%10.4g ", *mat00);
  printf("%10.4g ", *mat01);
  printf("%10.4g\n", *mat02);

  printf("%10.4g ", *mat10);
  printf("%10.4g ", *mat11);
  printf("%10.4g\n", *mat12);

  printf("%10.4g ", *mat20);
  printf("%10.4g ", *mat21);
  printf("%10.4g\n", *mat22);

}

/** print a vector
 * \param[in] double* v
 */
void print3(double* v)
{
  printf("%10.4g\n", *v++);
  printf("%10.4g\n", *v++);
  printf("%10.4g\n", *v);
}
