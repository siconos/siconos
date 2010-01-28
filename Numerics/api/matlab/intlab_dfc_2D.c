/*!\file intlab_dfc_2D.c


  This subroutine allows the dual resolution of contact problems with friction\n

  These problems are solved thanks to @ref dfc_2D solvers or @ref lcp solvers via a reformulation or a condensation.

  This subroutine is a matlab interface that allows you to use an @ref dfc_2D solver or @ref lcp solvers.

  You only have to enter


  \param K1           The stiffness, a vector of double (in which the components
                       of the matrix have a Fortran storage),

  \param F1           The right hand side, a vector of double,

  \param n            The dimension of the DFC_2D problem, an integer,

  \param mu           The friction coefficient, a double nonnegative,


  and a matlab structure with this following arguments:


  \param name         The name of the subroutine you want to use,

  \param itermax      The maximum number of iteration required , an integer,

  \param tol          The tolerance required, a double,

  \param k_latin      The research parameter, a double stricttly positive,

  \param chat         The chattering, an integer ( = 0: no comments, > 0 : comments)

  \param J1       The gap in normal contact direction, vector of doubles,

  \param dim_t        The dimension of ddl_t (= dimension of ddl_n), integer,

  \param ddl_n        The contact in normal direction dof (not prescribed), vector of integers,

  \param ddl_t        The contact in tangential direction dof (not prescribed),vector of integers,

  \param dim_d        The dimension of ddl_d, an integer

  \param ddl_d        The prescribed dof, vector of integers.


See go_dfc_2D.m for sample.


\author Nineb Sheherazade.
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"
#include "SiconosNumerics_lab.h"
#include "../../src/NSSpack/dfc_2D_solver.c"
#include "../../src/NSSpack/dfc_2D2cond_2D.c"
#include "../../src/NSSpack/cond_2D2dfc_2D.c"
#include "../../src/NSSpack/dfc_2D_latin.c"
#include "../../src/NSSpack/dfc_2D2lcp.c"
#include "../../src/NSSpack/lcp2dfc_2D.c"
#include "../../src/NSSpack/lcp_nlgs.c"
#include "../../src/NSSpack/lcp_latin.c"
#include "../../src/NSSpack/lcp_lexicolemke.c"
#include "../../src/NSSpack/lcp_lemke.c"
#include "../../src/NSSpack/diffns.c"





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
  double            *sort2, sort3, *J1;
  int               sort3bis, *pr3bis;
  double            *sort4;
  double            *sort11, *sort22, *sort33;
  int               *sort55, *pr5, *ddl_n, *ddl_tt, *ddl_d;


  double            tol, k_lat;
  int               itermax, chat ;

  char              *mot1, *mot2, *mot3, *mot4, *mot5, *mot6, *mot7;
  method            meth_dfc_2D;
  double            *vec, *q, *z, *w, *info, *mu;
  int               n, *nr, *dim_tt, *dim_d;



  mot1 = "NLGS";
  mot2 = "Cfd_latin";
  mot3 = "Lemke";
  /*  mot4 = "CPG";
      mot5 = "Latin";
      mot6 = "QP";

  */
  mot7 = "NSQP";



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



  meth_dfc_2D.dfc_2D.mu = *mu;



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

    meth_dfc_2D.dfc_2D.name = buf;

    /*                        2nd element of the structure                     */

    if ((strcmp(buf, mot2) == 0))
    {

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





      total_num_of_elements2 = 1;


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {
        if (field_index == 1)
        {

          itermax = pr3bis[index2];

          meth_dfc_2D.dfc_2D.itermax =  itermax;





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
          meth_dfc_2D.dfc_2D.tol = tol;



        }
      }


      /*                       End 3  element of the structure                     */
      /*                       4  element of the structure                     */

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

      if (category != mxDOUBLE_CLASS)
      {
        mexPrintf("The 4 element of the structure must be a DOUBLE");
      }


      if (field_index == 3)  pr3 = mxGetPr(field_array_ptr);


      total_num_of_elements2 = 1; /* mxGetNumberOfElements(field_array_ptr);*/


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      sort22 = (double*) malloc(1 * 1 * sizeof(double));



      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 3)
        {
          sort22[index2] = pr3[index2];

          k_lat = sort22[index2];
          meth_dfc_2D.dfc_2D.k_latin = k_lat;





        }
      }


      /*                      End of 4  element of the structure                     */

      /*                       5  element of the structure                     */



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


      if (category != mxINT32_CLASS)
      {
        mexPrintf("The 5 element of the structure must be an integer");
      }




      pr5 = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = 1; /*mxGetNumberOfElements(field_array_ptr);*/


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      sort55 = (int*) malloc(1 * 1 * sizeof(int));


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 4)
        {
          sort55[index2] = pr5[index2];

          chat = sort55[index2];
          meth_dfc_2D.dfc_2D.chat = chat;




        }
      }


      /*                       End of 5 element of the structure                     */

      /*                       6  element of the structure                     */

      field_index = 5;

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
        mexPrintf("The 6 element of the structure must be a DOUBLE");
      }


      if (field_index == 5)  J1 = mxGetPr(field_array_ptr);


      total_num_of_elements2 = mrows1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      meth_dfc_2D.dfc_2D.J1 = (double *) malloc(mrows1 * sizeof(double));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 5)
        {
          meth_dfc_2D.dfc_2D.J1[index2] = J1[index2];


        }
      }


      /*                      End of 6  element of the structure                     */

      /*                       7  element of the structure                     */



      field_index = 6;
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
        mexPrintf("The 7 element of the structure must be an integer");
      }




      dim_tt = mxGetData(field_array_ptr);


      total_num_of_elements2 = 1; /*mxGetNumberOfElements(field_array_ptr);*/


      meth_dfc_2D.dfc_2D.dim_tt = *dim_tt;/* (int *) malloc(*dim_tt*sizeof(int));*/


      /*                       End of 7 element of the structure                     */
      /*                       8  element of the structure                     */



      field_index = 7;
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
        mexPrintf("The 8 element of the structure must be an integer");
      }




      ddl_n = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_tt; /*mxGetNumberOfElements(field_array_ptr);*/





      meth_dfc_2D.dfc_2D.ddl_n = (int *) malloc(*dim_tt * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 7)
        {
          meth_dfc_2D.dfc_2D.ddl_n[index2] = ddl_n[index2] - 1;

        }
      }


      /*                       End of 8 element of the structure                     */

      /*                       9  element of the structure                     */



      field_index = 8;
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
        mexPrintf("The 9 element of the structure must be an integer");
      }




      ddl_tt = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_tt; /*mxGetNumberOfElements(field_array_ptr);*/



      meth_dfc_2D.dfc_2D.ddl_tt = (int *) malloc(*dim_tt * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 8)
        {
          meth_dfc_2D.dfc_2D.ddl_tt[index2] = ddl_tt[index2] - 1;



        }
      }


      /*                       End of 9 element of the structure                     */

      /*                       10  element of the structure                     */



      field_index = 9;
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
        mexPrintf("The 10 element of the structure must be an integer");
      }




      dim_d = mxGetData(field_array_ptr);


      total_num_of_elements2 = 1; /*mxGetNumberOfElements(field_array_ptr);*/




      meth_dfc_2D.dfc_2D.dim_d = *dim_d;/* (int *) malloc(*dim_tt*sizeof(int));*/



      /*                       End of 10 element of the structure                     */
      /*                       11  element of the structure                     */



      field_index = 10;
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
        mexPrintf("The 11 element of the structure must be an integer");
      }




      ddl_d = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_d;






      meth_dfc_2D.dfc_2D.ddl_d = (int *) malloc(*dim_d * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 10)
        {
          meth_dfc_2D.dfc_2D.ddl_d[index2] = ddl_d[index2] - 1;






        }
      }

      free(sort2);
      free(sort22);
      free(sort55);
    }

    /*                       End of 11 element of the structure                     */




    if ((strcmp(buf, mot1) == 0))
    {

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





      total_num_of_elements2 = 1;


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {
        if (field_index == 1)
        {

          itermax = pr3bis[index2];

          meth_dfc_2D.dfc_2D.itermax =  itermax;





        }
      }


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
          meth_dfc_2D.dfc_2D.tol = tol;



        }
      }


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
          meth_dfc_2D.dfc_2D.chat = chat;




        }
      }


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


      if (field_index == 4)  J1 = mxGetPr(field_array_ptr);


      total_num_of_elements2 = mrows1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      meth_dfc_2D.dfc_2D.J1 = (double *) malloc(mrows1 * sizeof(double));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 4)
        {
          meth_dfc_2D.dfc_2D.J1[index2] = J1[index2];


        }
      }



      field_index = 5;
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
        mexPrintf("The 6 element of the structure must be an integer");
      }




      dim_tt = mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;


      meth_dfc_2D.dfc_2D.dim_tt = *dim_tt;



      field_index = 6;
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
        mexPrintf("The 7 element of the structure must be an integer");
      }




      ddl_n = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_tt;





      meth_dfc_2D.dfc_2D.ddl_n = (int *) malloc(*dim_tt * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 6)
        {
          meth_dfc_2D.dfc_2D.ddl_n[index2] = ddl_n[index2] - 1;

        }
      }



      field_index = 7;
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
        mexPrintf("The 8 element of the structure must be an integer");
      }




      ddl_tt = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_tt;



      meth_dfc_2D.dfc_2D.ddl_tt = (int *) malloc(*dim_tt * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 7)
        {
          meth_dfc_2D.dfc_2D.ddl_tt[index2] = ddl_tt[index2] - 1;



        }
      }


      field_index = 8;
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
        mexPrintf("The 9 element of the structure must be an integer");
      }




      dim_d = mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;




      meth_dfc_2D.dfc_2D.dim_d = *dim_d;


      field_index = 9;
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
        mexPrintf("The 10 element of the structure must be an integer");
      }




      ddl_d = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_d;






      meth_dfc_2D.dfc_2D.ddl_d = (int *) malloc(*dim_d * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 9)
        {
          meth_dfc_2D.dfc_2D.ddl_d[index2] = ddl_d[index2] - 1;
        }
      }

      free(sort2);
      free(sort55);
    }



    if ((strcmp(buf, mot3) == 0))
    {

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





      total_num_of_elements2 = 1;


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {
        if (field_index == 1)
        {

          itermax = pr3bis[index2];

          meth_dfc_2D.dfc_2D.itermax =  itermax;
        }
      }


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


      if (category != mxINT32_CLASS)
      {
        mexPrintf("The 3 element of the structure must be an integer");
      }




      pr5 = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      sort55 = (int*) malloc(1 * 1 * sizeof(int));


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 2)
        {
          sort55[index2] = pr5[index2];

          chat = sort55[index2];
          meth_dfc_2D.dfc_2D.chat = chat;




        }
      }


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

      if (category != mxDOUBLE_CLASS)
      {
        mexPrintf("The 4 element of the structure must be a DOUBLE");
      }


      if (field_index == 3)  J1 = mxGetPr(field_array_ptr);


      total_num_of_elements2 = mrows1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      meth_dfc_2D.dfc_2D.J1 = (double *) malloc(mrows1 * sizeof(double));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 3)
        {
          meth_dfc_2D.dfc_2D.J1[index2] = J1[index2];


        }
      }



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


      if (category != mxINT32_CLASS)
      {
        mexPrintf("The 5 element of the structure must be an integer");
      }




      dim_tt = mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;


      meth_dfc_2D.dfc_2D.dim_tt = *dim_tt;



      field_index = 5;
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
        mexPrintf("The 6 element of the structure must be an integer");
      }




      ddl_n = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_tt;





      meth_dfc_2D.dfc_2D.ddl_n = (int *) malloc(*dim_tt * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 5)
        {
          meth_dfc_2D.dfc_2D.ddl_n[index2] = ddl_n[index2] - 1;

        }
      }



      field_index = 6;
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
        mexPrintf("The 7 element of the structure must be an integer");
      }




      ddl_tt = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_tt;



      meth_dfc_2D.dfc_2D.ddl_tt = (int *) malloc(*dim_tt * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 6)
        {
          meth_dfc_2D.dfc_2D.ddl_tt[index2] = ddl_tt[index2] - 1;



        }
      }


      field_index = 7;
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
        mexPrintf("The 8 element of the structure must be an integer");
      }




      dim_d = mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;




      meth_dfc_2D.dfc_2D.dim_d = *dim_d;


      field_index = 8;
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
        mexPrintf("The 9 element of the structure must be an integer");
      }




      ddl_d = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_d;






      meth_dfc_2D.dfc_2D.ddl_d = (int *) malloc(*dim_d * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 8)
        {
          meth_dfc_2D.dfc_2D.ddl_d[index2] = ddl_d[index2] - 1;
        }
      }

      free(sort55);
    }


    if ((strcmp(buf, mot7) == 0))
    {

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

      if (category != mxDOUBLE_CLASS)
      {
        mexPrintf("The 2 element of the structure must be a DOUBLE");
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
          meth_dfc_2D.dfc_2D.tol = tol;



        }
      }


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


      if (category != mxINT32_CLASS)
      {
        mexPrintf("The 3 element of the structure must be an integer");
      }




      pr5 = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      sort55 = (int*) malloc(1 * 1 * sizeof(int));


      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 2)
        {
          sort55[index2] = pr5[index2];

          chat = sort55[index2];
          meth_dfc_2D.dfc_2D.chat = chat;




        }
      }


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

      if (category != mxDOUBLE_CLASS)
      {
        mexPrintf("The 4 element of the structure must be a DOUBLE");
      }


      if (field_index == 3)  J1 = mxGetPr(field_array_ptr);


      total_num_of_elements2 = mrows1;


      mrows = mxGetM(field_array_ptr);
      mcols = mxGetN(field_array_ptr);


      meth_dfc_2D.dfc_2D.J1 = (double *) malloc(mrows1 * sizeof(double));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 3)
        {
          meth_dfc_2D.dfc_2D.J1[index2] = J1[index2];


        }
      }



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


      if (category != mxINT32_CLASS)
      {
        mexPrintf("The 5 element of the structure must be an integer");
      }




      dim_tt = mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;


      meth_dfc_2D.dfc_2D.dim_tt = *dim_tt;



      field_index = 5;
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
        mexPrintf("The 6 element of the structure must be an integer");
      }




      ddl_n = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_tt;





      meth_dfc_2D.dfc_2D.ddl_n = (int *) malloc(*dim_tt * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 5)
        {
          meth_dfc_2D.dfc_2D.ddl_n[index2] = ddl_n[index2] - 1;

        }
      }



      field_index = 6;
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
        mexPrintf("The 7 element of the structure must be an integer");
      }




      ddl_tt = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_tt;



      meth_dfc_2D.dfc_2D.ddl_tt = (int *) malloc(*dim_tt * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 6)
        {
          meth_dfc_2D.dfc_2D.ddl_tt[index2] = ddl_tt[index2] - 1;



        }
      }


      field_index = 7;
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
        mexPrintf("The 8 element of the structure must be an integer");
      }




      dim_d = mxGetData(field_array_ptr);


      total_num_of_elements2 = 1;




      meth_dfc_2D.dfc_2D.dim_d = *dim_d;


      field_index = 8;
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
        mexPrintf("The 9 element of the structure must be an integer");
      }




      ddl_d = (int *) mxGetData(field_array_ptr);


      total_num_of_elements2 = *dim_d;






      meth_dfc_2D.dfc_2D.ddl_d = (int *) malloc(*dim_d * sizeof(int));

      for (index2 = 0; index2 < total_num_of_elements2; index2++)
      {

        if (field_index == 8)
        {
          meth_dfc_2D.dfc_2D.ddl_d[index2] = ddl_d[index2] - 1;
        }
      }

      free(sort55);
    }



    printf("we enter in dfc_2D_solver using the %s method \n\n", buf);

    if (strcmp(buf, mot1) == 0)
    {
      toto = dfc_2D_solver(vec, q, &mrows1, &meth_dfc_2D, z, w);
    }
    else if (strcmp(buf, mot2) == 0)
      toto = dfc_2D_solver(vec, q, &mrows1, &meth_dfc_2D, z, w);
    else if (strcmp(buf, mot3) == 0)
      toto = dfc_2D_solver(vec, q, &mrows1, &meth_dfc_2D, z, w);
    else if (strcmp(buf, mot7) == 0)
      toto = dfc_2D_solver(vec, q, &mrows1, &meth_dfc_2D, z, w);
    /*   else if (strcmp(buf, mot5) == 0)
     toto = dfc_2D_solver(vec, q, &mrows1, &meth_dfc_2D, z, w);*/
    else printf("Warning : Unknown solving method : %s\n", buf);

    *info = toto;





    free(meth_dfc_2D.dfc_2D.J1);
    free(meth_dfc_2D.dfc_2D.ddl_n);
    free(meth_dfc_2D.dfc_2D.ddl_tt);
    free(meth_dfc_2D.dfc_2D.ddl_d);


  }
}
