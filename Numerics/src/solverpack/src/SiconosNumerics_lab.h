#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
  char * nom_method;
  int itermax;
  double tol;
  double mu;
  double k_latin;
} methode_cfp;


typedef struct
{
  char *name;
  int itermax;
  double tol;
  double mu;
  double k_latin;
  double *J1;
  int *ddl_i;
  int *ddl_n;
  int *ddl_tt;
  int *ddl_c;
  int dim_i;
  int dim_c;
  int dim_tt;
  int dim_n;
} method_dfc_2D;


typedef struct
{
  char * name;
  int itermax;
  double tol;
  double k_latin;
} method_lcp;


typedef struct
{
  char * nom_method;
  int itermax;
  double tol;
  double *a;
  double *b;
  double k_latin;
} methode_rp;



typedef struct
{
  char * nom_method;
  int itermax;
  double tol;
  double *a;
  double *b;
  double k_latin;
} methode_rd;




typedef union
{
  methode_cfp cfp;
  methode_cfd cfd;
  methode_lcp lcp;
  methode_rp rp;
  methode_rd rd;

} methode;





