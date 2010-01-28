#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct
{

  char    * name;
  int     itermax;
  int     chat;
  double  tol;
  double  mu;
  double  k_latin;
  char    *normType;
  int     iter;
  double  err;

} method_pfc_2D;


typedef struct
{

  char    *name;
  int     itermax;
  double  tol;
  double  mu;
  double  k_latin;
  double  *J1;
  int     *ddl_n;
  int     *ddl_tt;
  int     *ddl_d;
  int     dim_tt;
  int     dim_d;
  int     chat;
  char    *normType;
  int     iter;
  double  err;

} method_dfc_2D;



typedef struct
{


  char    *name;
  int     itermax;
  double  tol;
  double  k_latin;
  double  relax;
  int     chat;
  char    *normType;
  int     iter;
  double  err;

} method_lcp;



typedef struct
{

  char    * name;
  int     itermax;
  double  tol;
  double  *a;
  double  *b;
  double  k_latin;
  int     chat;
  char    *normType;
  int     iter;
  double  err;

} method_pr;



typedef struct
{

  char    *name;
  int     itermax;
  double  tol;
  double  k_latin;
  double  *a;
  double  *b;
  int     chat;
  char    *normType;
  int     iter;
  double  err;

} method_dr;




typedef union
{

  method_pfc_2D  pfc_2D;
  method_dfc_2D  dfc_2D;
  method_lcp     lcp;
  method_pr      pr;
  method_dr      dr;

} method;





/*int lcp_solver( double *vec, double *q , int *n , method *pt , double *z , double *w );


void lcp_lemke( double *vec, double *qqq,int *nn, int *itermax, double *zlem,
      double *wlem, int *it_end, double *res, int *info );

void lcp_qp( int *nn , double *vec , double *q , double *z , double *w , int *info ,
    int *iparamLCP , double *dparamLCP );

void lcp_cpg( int *nn , double *vec , double *q , double *z , double *w , int *info ,
     int *iparamLCP , double *dparamLCP );

void lcp_nlgs( int *nn , double *vec , double *q , double *z , double *w , int *info ,
      int *iparamLCP , double *dparamLCP );

void lcp_nsqp( int *nn , double *vec , double *q , double *z , double *w , int *info ,
      int *iparamLCP , double *dparamLCP );

void lcp_latin( int *nn , double *vec , double *q , double *z , double *w , int *info ,
       int *iparamLCP , double *dparamLCP );


void lcp_latin_w( int *nn , double *vec , double *q , double *z , double *w , int *info ,
       int *iparamLCP , double *dparamLCP );


void lcp_lexicolemke( int *nn , double *vec , double *q , double *z , double *w , int *info ,
      int *iparamLCP , double *dparamLCP );

void lcp_newton_min( int *nn , double *vec , double *q , double *z , double *w , int *info ,
      int *iparamLCP , double *dparamLCP );


*/
