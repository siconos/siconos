/*!\file intlab_lcp.c


   This subroutine allows the resolution of LCP (Linear Complementary Problem).
   Try \f$(z,w)\f$ such that:

\f$
\left\lbrace
\begin{array}{l}
w - M z = q\\
0 \le z \perp w \ge 0\\
\end{array}
\right.
\f$

  here M is an n by n  matrix, q an n-dimensional vector, w an n-dimensional  vector and z an n-dimensional vector.
  This system of equalities and inequalities is solved thanks to @ref lcp solvers.

  This subroutine is a matlab interface that allows you to use an @ref lcp solver.

  You only have to enter the M matrix, the q vector, the dimension row of q, and a matlab structure with this following arguments:

  * the name of the solver 'NLGS', 'LexicoLemke', 'CPG', 'Latin','Latin_w', 'QP', 'NSQP' or 'NewtonMin'.
  * the maximum of iterations required,
  * the tolerance value,
  * the search direction for the Latin,
  * the relaxation parameter,
  * the chattering.


See go_lcp.m for sample.

\author Nineb Sheherazade.
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include "SiconosNumerics_lab.h"
#include "../../src/NSSpack/lcp_solver.c"
#include "../../src/NSSpack/lcp_nlgs.c"
#include "../../src/NSSpack/lcp_latin.c"
#include "../../src/NSSpack/lcp_lemke.c"
#include "../../src/NSSpack/lcp_cpg.c"
#include "../../src/NSSpack/lcp_lexicolemke.c"
#include "../../src/NSSpack/lcp_qp.c"
#include "../../src/NSSpack/lcp_latin_w.c"
#include "../../src/NSSpack/lcp_nsqp.c"
#include "../../src/NSSpack/lcp_newton_min.c"



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
  double            *sort11, *sort22, *sort33, *sort44;
  int               *sort55, *pr5, *pr33bis, *sort22bis, *pr22bis;


  double            tol, k_lat, relax;
  int               itermax, chat;

  char              *mot1, *mot2, *mot3, *mot4, *mot5, *mot6, *mot7, *mot8;
  method            meth_lcp;
  double            *vec, *q, *z, *w, *info, *mu;
  int               n, *nr;



  mot1 = "NLGS";
  mot2 = "Latin";
  mot3 = "CPG";
  mot4 = "Latin_w";
  mot5 = "LexicoLemke";
  mot6 = "NewtonMin";
  mot7 = "QP";
  mot8 = "NSQP";







  if (nlhs != 3)
  {
    mexErrMsgTxt("3 output required.");
  }


  if (nrhs != 4)
  {
    mexErrMsgTxt("4 input required.");
  }


  vec = mxGetPr(prhs[0]);
  q   = mxGetPr(prhs[1]);
  nr  = mxGetData(prhs[2]);



  n   = *nr;



  mrows1 = mxGetM(prhs[1]);
  mcols1 = mxGetN(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(mrows1, mcols1, mxREAL); /* z */
  plhs[1] = mxCreateDoubleMatrix(mrows1, mcols1, mxREAL); /* w */
  plhs[2] = mxCreateDoubleMatrix(mcols1, mcols1, mxREAL); /* info */


  z = mxGetPr(plhs[0]);
  w = mxGetPr(plhs[1]);
  info = mxGetPr(plhs[2]);

  category = mxGetClassID(prhs[3]);

  if (category != mxSTRUCT_CLASS)
  {
    mexPrintf("The thrid input must be a structure");
  }


  total_num_of_elements = mxGetNumberOfElements(prhs[3]);
  number_of_fields = mxGetNumberOfFields(prhs[3]);


  mexPrintf("\n\t\t");

  index = 0;
  field_index = 0;

  field_array_ptr = mxGetFieldByNumber(prhs[3],
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

    meth_lcp.lcp.name = buf;




    if ((strcmp(buf, mot1) == 0) || (strcmp(buf, mot2) == 0) || (strcmp(buf, mot3) == 0) || (strcmp(buf, mot4) == 0) || (strcmp(buf, mot6) == 0))

    {

      /*                        2nd element of the structure                     */

      if (number_of_fields < 4)
        mexPrintf("Number of elements in the structure not valid.");


      field_index = 1;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
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





      total_num_of_elements2 = 1;


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {
        if (field_index == 1)
        {
          itermax = pr3bis[index2];
          meth_lcp.lcp.itermax =  itermax;

        }
      }

      /*                      End of  2nd element of the structure                     */

      /*                        3 element of the structure                     */


      field_index = 2;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
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


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      sort2 = (double*) malloc(1 * 1 * sizeof(double));



      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 2)
        {


          sort2[index2] = pr2[index2];

          tol = sort2[index2];
          meth_lcp.lcp.tol = tol;



        }
      }

      /*                       End 3  element of the structure                     */


      /*                       4 element of the structure                     */



      field_index = 3;
      field_array_ptr = mxGetFieldByNumber(prhs[3],
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


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);

      sort55 = (int*) malloc(1 * 1 * sizeof(int));


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 3)
        {
          sort55[index2] = pr5[index2];

          chat = sort55[index2];
          meth_lcp.lcp.chat = chat;




        }
      }

      free(sort2);
      free(sort55);
    }



    if (strcmp(buf, mot2) == 0)
    {

      number_of_fields = mxGetNumberOfFields(prhs[3]);


      if (number_of_fields != 5)
        mexPrintf("Number of elements in the structure not valid.");




      field_index = 4;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
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


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      sort22 = (double*) malloc(1 * 1 * sizeof(double));



      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 4)
        {
          sort22[index2] = pr3[index2];

          k_lat = sort22[index2];
          meth_lcp.lcp.k_latin = k_lat;





        }
      }

      free(sort22);

    }


    if (strcmp(buf, mot4) == 0)
    {


      number_of_fields = mxGetNumberOfFields(prhs[3]);


      if (number_of_fields != 6)
        mexPrintf("Number of elements in the structure not valid.");



      field_index = 4;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
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


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);

      sort22 = (double*) malloc(1 * 1 * sizeof(double));



      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 4)
        {
          sort22[index2] = pr3[index2];

          k_lat = sort22[index2];
          meth_lcp.lcp.k_latin = k_lat;

        }
      }

      free(sort22);




      field_index = 5;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
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
        mexPrintf("The 6 element of the structure must be a DOUBLE");
      }


      if (field_index == 5)  pr4 = mxGetPr(field_array_ptr);


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);

      sort44 = (double*) malloc(1 * 1 * sizeof(double));



      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 5)
        {
          sort44[index2] = pr4[index2];

          relax = sort44[index2];

          meth_lcp.lcp.relax = relax;


        }
      }


      free(sort44);

    }



    if ((strcmp(buf, mot5) == 0))

    {

      field_index = 1;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
                                           index,
                                           field_index);
      if (number_of_fields < 3)
        mexPrintf("Number of elements in the structure not valid.");


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






      pr33bis = (int *) mxGetData(field_array_ptr);





      total_num_of_elements2 = 1;


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {
        if (field_index == 1)
        {

          itermax = pr33bis[index2];
          meth_lcp.lcp.itermax =  itermax;





        }
      }

      field_index = 2;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
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
        mexPrintf("The 2 element of the structure must be an INTEGER");
      }

      if (field_index == 2) pr22bis = mxGetPr(field_array_ptr);


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);

      sort22bis = (int*) malloc(1 * 1 * sizeof(double));



      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 2)
        {


          sort22bis[index2] = pr22bis[index2];

          chat = sort22bis[index2];
          meth_lcp.lcp.chat = chat;



        }
      }

      free(sort22bis);


    }

    if ((strcmp(buf, mot7) == 0) || (strcmp(buf, mot8) == 0))

    {

      field_index = 1;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
                                           index,
                                           field_index);
      if (number_of_fields < 3)
        mexPrintf("Number of elements in the structure not valid.");


      if (field_array_ptr == NULL)
        mexPrintf("\tEmpty Field\n");
      else
      {
        category = mxGetClassID(field_array_ptr);
      }




      if (category != mxDOUBLE_CLASS)
      {
        mexPrintf("The 2 element of the structure must be a double\n");
      }


      pr2 = (double *) mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {
        if (field_index == 1)
        {

          tol = pr2[index2];
          meth_lcp.lcp.tol = tol;


        }
      }


      field_index = 2;

      field_array_ptr = mxGetFieldByNumber(prhs[3],
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
        mexPrintf("The 2 element of the structure must be an INTEGER");
      }

      if (field_index == 2) pr22bis = mxGetPr(field_array_ptr);


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);

      sort22bis = (int*) malloc(1 * 1 * sizeof(double));



      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 2)
        {


          sort22bis[index2] = pr22bis[index2];

          chat = sort22bis[index2];
          meth_lcp.lcp.chat = chat;



        }
      }

      free(sort22bis);


    }


    printf("we enter in lpc_solver using the %s method \n\n", buf);


    if (strcmp(buf, mot1) == 0)
    {
      toto = lcp_solver(vec, q, &mrows1, &meth_lcp, z, w);
    }
    else if (strcmp(buf, mot2) == 0)
      toto = lcp_solver(vec, q, &mrows1, &meth_lcp, z, w);
    else if (strcmp(buf, mot3) == 0)
      toto = lcp_solver(vec, q, &mrows1, &meth_lcp, z, w);
    else if (strcmp(buf, mot4) == 0)
      toto = lcp_solver(vec, q, &mrows1, &meth_lcp, z, w);
    else if (strcmp(buf, mot5) == 0)
      toto = lcp_solver(vec, q, &mrows1, &meth_lcp, z, w);
    else if (strcmp(buf, mot6) == 0)
      toto = lcp_solver(vec, q, &mrows1, &meth_lcp, z, w);
    else if (strcmp(buf, mot7) == 0)
      toto = lcp_solver(vec, q, &mrows1, &meth_lcp, z, w);
    else if (strcmp(buf, mot8) == 0)
      toto = lcp_solver(vec, q, &mrows1, &meth_lcp, z, w);
    else printf("Warning : Unknown solving method : %s\n", buf);

    *info = toto;





  }
}
