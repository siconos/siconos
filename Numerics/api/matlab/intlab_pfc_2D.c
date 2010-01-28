/*!\file intlab_pfc_2D.c

 This subroutine allows the primal resolution of contact problems with friction.

Try \f$(z,w)\f$ such that:
\f$
\left\lbrace
\begin{array}{l}
w - M z=q\\
0 \le z_n \perp w_n \ge 0\\
-w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
\end{array}
\right.
\f$

 here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.

 This system of equations and inequalities is solved thanks to @ref pfc_2D solvers.
 The routine's call is due to the function pfc_2D_solver.c.

  This subroutine is a matlab interface that allows you to use an @ref pfc_2D solver.

  You only have to enter the M matrix, the q vector, the dimension row of q, the friction coefficient mu and a matlab structure with this following arguments:

  * the name of the solver 'NLGS', 'CPG' or 'Latin',
  * the maximum of iterations required,
  * the tolerance value,
  * the search direction for the Latin,
  * the chattering.

See go_pfc_2D.m for sample.


\author Nineb Sheherazade.
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include "SiconosNumerics_lab.h"
#include "../../src/NSSpack/pfc_2D_solver.c"
#include "../../src/NSSpack/pfc_2D_latin.c"
#include "../../src/NSSpack/pfc_2D_cpg.c"
#include "../../src/NSSpack/pfc_2D_nlgs.c"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int               i;
  mxClassID         category;
  int               total_num_of_elements, number_of_fields, index, field_index;
  const char        *field_name;

  const mxArray     *field_array_ptr;
  char              *buf, *sort1;
  int               number_of_dimensions;
  const int         *dims;
  int               buflen, d, page, total_number_of_pages, elements_per_page;
  double            *pr, *pr2, *pr3, *pr4;
  int               total_num_of_elements2, index2;

  int               mrows, mcols, mrows1, mcols1, toto;
  double            *sort2, sort3;
  int               sort3bis, *pr3bis;
  double            *sort4;
  double            *sort11, *sort22, *sort33;
  int               *sort55, *pr5;


  double            tol, k_lat;
  int               itermax, chat ;

  char              *mot1, *mot2, *mot3;
  method            meth_pfc_2D;
  double            *vec, *q, *z, *w, *info, *mu;
  int               n, *nr;



  mot1 = "NLGS";
  mot2 = "Latin";
  mot3 = "CPG";




  if (nlhs != 3)
  {
    mexErrMsgTxt("3 output required.");
  }


  if (nrhs != 5)
  {
    mexErrMsgTxt("5 input required.");
  }


  vec = mxGetPr(prhs[0]);
  q   = mxGetPr(prhs[1]);
  nr  = mxGetData(prhs[2]);
  mu  = mxGetPr(prhs[3]);

  meth_pfc_2D.pfc_2D.mu = *mu;



  n   = *nr;



  mrows1 = mxGetM(prhs[1]);
  mcols1 = mxGetN(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(mrows1, mcols1, mxREAL); /* z */
  plhs[1] = mxCreateDoubleMatrix(mrows1, mcols1, mxREAL); /* w */
  plhs[2] = mxCreateDoubleMatrix(mcols1, mcols1, mxREAL); /* info */


  z = mxGetPr(plhs[0]);
  w = mxGetPr(plhs[1]);
  info = mxGetPr(plhs[2]);

  category = mxGetClassID(prhs[4]);

  if (category != mxSTRUCT_CLASS)
  {
    mexPrintf("The thrid input must be a structure");
  }


  total_num_of_elements = mxGetNumberOfElements(prhs[4]);
  number_of_fields = mxGetNumberOfFields(prhs[4]);

  if (number_of_fields < 4)
    mexPrintf("Number of elements in the structure not valid.");



  mexPrintf("\n\t\t");

  index = 0;
  field_index = 0;

  field_array_ptr = mxGetFieldByNumber(prhs[4],
                                       index,
                                       field_index);


  if (field_array_ptr == NULL)
    mexPrintf("\tEmpty Field\n");
  else
  {


    category = mxGetClassID(field_array_ptr);

    if (category != mxCHAR_CLASS)
    {
      mexPrintf("The first element of the structure must be a CHAR");
    }



    buflen = mxGetNumberOfElements(field_array_ptr) + 1;
    buf = mxCalloc(buflen, sizeof(char));

    /* Copy the string data from string_array_ptr and place it into buf. */
    if (mxGetString(field_array_ptr, buf, buflen) != 0)
      mexErrMsgTxt("Could not convert string data.");

    printf("\n");

    meth_pfc_2D.pfc_2D.name = buf;

    /*                        2nd element of the structure                     */

    field_index = 1;

    field_array_ptr = mxGetFieldByNumber(prhs[4],
                                         index,
                                         field_index);


    if (field_array_ptr == NULL)
      mexPrintf("\tEmpty Field\n");
    else
    {
      category = mxGetClassID(field_array_ptr);
    }




    if (category != mxINT32_CLASS)
    {
      mexPrintf("The 2 element of the structure must be an integer\n");
    }






    pr3bis = (int *) mxGetData(field_array_ptr);





    total_num_of_elements2 = 1; /* mxGetNumberOfElements(field_array_ptr);*/


    for (index2 = 0; index2 < total_num_of_elements2; index2++)
    {
      if (field_index == 1)
      {

        /* sort3 = pr[index2];
           sort3 = pr[index2]; */

        itermax = pr3bis[index2];
        /*itermax = sort3;*/
        meth_pfc_2D.pfc_2D.itermax =  itermax;





      }
    }

    /*                      End of  2nd element of the structure                     */

    /*                        3 element of the structure                     */
    field_index = 2;

    field_array_ptr = mxGetFieldByNumber(prhs[4],
                                         index,
                                         field_index);


    if (field_array_ptr == NULL)
      mexPrintf("\tEmpty Field\n");
    else
    {
      category = mxGetClassID(field_array_ptr);
    }

    if (category != mxDOUBLE_CLASS)
    {
      mexPrintf("The 3 element of the structure must be a DOUBLE");
    }

    if (field_index == 2) pr2 = mxGetPr(field_array_ptr);


    total_num_of_elements2 = 1;/*mxGetNumberOfElements(field_array_ptr);*/


    mrows = mxGetM(field_array_ptr);
    mcols = mxGetN(field_array_ptr);


    /*      sort2 = (double*) malloc (mrows*mcols*sizeof(double));*/
    sort2 = (double*) malloc(1 * 1 * sizeof(double));



    for (index2 = 0; index2 < total_num_of_elements2; index2++)
    {

      if (field_index == 2)
      {


        sort2[index2] = pr2[index2];

        tol = sort2[index2];
        meth_pfc_2D.pfc_2D.tol = tol;



      }
    }


    /*                       End 3  element of the structure                     */


    /*                       5  element of the structure                     */



    field_index = 3;
    field_array_ptr = mxGetFieldByNumber(prhs[4],
                                         index,
                                         field_index);


    if (field_array_ptr == NULL)
      mexPrintf("\tEmpty Field\n");
    else
    {
      category = mxGetClassID(field_array_ptr);
    }


    if (category != mxINT32_CLASS)
    {
      mexPrintf("The 4 element of the structure must be an integer");
    }




    pr5 = (int *) mxGetData(field_array_ptr);


    total_num_of_elements2 = 1; /*mxGetNumberOfElements(field_array_ptr);*/


    mrows = mxGetM(field_array_ptr);
    mcols = mxGetN(field_array_ptr);


    /*      sort55 = (int*) malloc (mrows*mcols*sizeof(int));*/
    sort55 = (int*) malloc(1 * 1 * sizeof(int));


    for (index2 = 0; index2 < total_num_of_elements2; index2++)
    {

      if (field_index == 3)
      {
        sort55[index2] = pr5[index2];

        chat = sort55[index2];
        meth_pfc_2D.pfc_2D.chat = chat;




      }
    }


    /*                       End of 5 element of the structure                     */


    /*                       4  element of the structure                     */


    if (strcmp(buf, mot2) == 0)
    {

      field_index = 4;

      field_array_ptr = mxGetFieldByNumber(prhs[4],
                                           index,
                                           field_index);


      if (field_array_ptr == NULL)
        mexPrintf("\tEmpty Field\n");
      else
      {
        category = mxGetClassID(field_array_ptr);
      }

      if (category != mxDOUBLE_CLASS)
      {
        mexPrintf("The 5 element of the structure must be a DOUBLE");
      }


      if (field_index == 4)  pr3 = mxGetPr(field_array_ptr);


      total_num_of_elements2 = 1; /* mxGetNumberOfElements(field_array_ptr);*/


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      /*      sort22 = (double*) malloc (mrows*mcols*sizeof(double));*/
      sort22 = (double*) malloc(1 * 1 * sizeof(double));



      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 4)
        {
          sort22[index2] = pr3[index2];

          k_lat = sort22[index2];
          meth_pfc_2D.pfc_2D.k_latin = k_lat;





        }
      }

      free(sort22);

    }
    /*                      End of 4  element of the structure                     */




    printf("we enter in pfc_2D_solver using the %s method \n\n", buf);


    if (strcmp(buf, mot1) == 0)
    {
      toto = pfc_2D_solver(vec, q, &mrows1, &meth_pfc_2D, z, w);
    }
    else if (strcmp(buf, mot2) == 0)
      toto = pfc_2D_solver(vec, q, &mrows1, &meth_pfc_2D, z, w);
    else if (strcmp(buf, mot3) == 0)
      toto = pfc_2D_solver(vec, q, &mrows1, &meth_pfc_2D, z, w);
    else printf("Warning : Unknown solving method : %s\n", buf);

    *info = toto;


    free(sort2);
    free(sort55);




  }
}
