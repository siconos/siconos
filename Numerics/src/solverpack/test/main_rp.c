///////////////////////////////////////////////////////////////////////////
//  This main file allows the primal resolution of contact problems with friction:
//  try (z,w) such that:
//
//
//
//                    M z  - w = q                 (1)
//                    0<= zn , 0<= wn , zn.wn= O   (2) (Signorini contact law)
//                    -wt in diff(PSI)   (zt)      (3) (Coulomb friction law)
//                                 [-mu*zn, mu*zn]
//
//
//  here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional//  vector and w an n-dimensional vector.
//
//  This system of equations and inequalities is solved thanks to cfp_subroutine:
//        cfp_gsnl (M,q,n,mu,itermax,tol,z,w,it_end,res,info)
//        cfp_gcp (M,q,n,mu,itermax,tol,z,w,it_end,res,info)
//        cfp_latin (M,q,n,k_latin,mu,itermax,tol,z,w,it_end,res,info)
//
//  where _ itermax is the maximum iterations required, it's an integer
//        _ mu is the friction coefficient, it's a float
//        _ res is the residue, it's a float
//        _ it_end is the number of iterations carried out, it's an integer
//        _ tol is the tolerance value desired, it's a float
//        _ k_latin is the parameter research of the latin, it's a float
//        _ z and w are the solutions of the problem
//        _ info shows the termination reason,0 is successful otherwise 1, it's an integer.
//
//
//    For more information about the methods see the chapter 4 of the Siconos manual theory.
//
//
//
//  The subroutine's call is due to the function solve_cfp:
//
//  int solve_cfp (double (*M)[maxcols],double *q,int n,methode *pt, double *z,double *w)
//
//  where M is an n by n matrix, q an n-dimensional vector, n is the row dimension
//  of M, and pt a pointer other a structure ( methode). z and w are n-dimensional
//  vectors solution.
//  methode is a variable with a structure type; this structure gives to the function
// solve_pcf, the name and the parameters (itermax, tol, k_latin,...) of the method we want to use.
//  This function return an interger:  0 successful return otherwise 1.
//
//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SiconosNumerics.h"


void test1(void)
{
  FILE *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9;
  int i, j, nl, nc, nll, it, info, n = 40, dimM = n, incx, incy;
  double *q, *z, *w, *vec, *a, *b, *c, *qqt, *zt;
  double(*M)[n], (*Mt)[n];
  double qi, Mij, alpha, beta;
  char val[14], vall[14];
  methode meth_rp;
  char trans;

  //  meth_rp.rp.nom_method="latin";
  strcpy(meth_rp.rp.nom_method, "Gsnl");
  meth_rp.rp.itermax = 50;
  meth_rp.rp.tol = 0.001;
  meth_rp.rp.k_latin = 0.003; //0.00005;



  if ((f1 = fopen("M_relay.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }



  M = malloc(dimM * dimM * sizeof(double));
  Mt = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));


  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////

    *(*(M + nc - 1) + nl - 1) = Mij;
    *(*(Mt + nl - 1) + nc - 1) = Mij;

    //////////////         fin transpos         ////////////////////

  }

  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  for (i = 0; i < dimM; i++)
  {
    for (j = 0; j < dimM; j++)
    {
      vec[j * dimM + i] = M[i][j];

    }
  }
  ////////////////////////////////////////////////////////////////////////




  if ((f2 = fopen("q_relay.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  if ((f5 = fopen("a_relay.dat", "r")) == NULL)
  {
    perror("fopen 5");
    exit(5);
  }


  if ((f6 = fopen("b_relay.dat", "r")) == NULL)
  {
    perror("fopen 6");
    exit(6);
  }



  q = malloc(dimM * sizeof(double));
  qqt = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  zt = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));
  a = malloc(dimM * sizeof(double));
  b = malloc(dimM * sizeof(double));
  c = malloc(dimM * sizeof(double));




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    *(q + nll - 1) = qi;
  }


  while (!feof(f5))
  {
    fscanf(f5, "%d", &nll);
    fscanf(f5, "%s", vall);
    qi = atof(vall);
    fprintf(f5, "%d %.14e\n", nll, qi);
    *(a + nll - 1) = qi;
  }


  while (!feof(f6))
  {
    fscanf(f6, "%d", &nll);
    fscanf(f6, "%s", vall);
    qi = atof(vall);
    fprintf(f6, "%d %.14e\n", nll, qi);
    *(b + nll - 1) = qi;
  }


  for (i = 0; i < dimM; i++)
  {
    c[i] = (-b[i] - a[i]) / 2;
    qqt[i] = q[i];
  }


  trans = 'T';
  alpha = 1.;
  beta = 1.;
  incx = 1;
  incy = 1;


  dgemv_(&trans, &n, &n, &alpha, M, &n, c, &incx, &beta, qqt, &incy);

  meth_rp.rp.a = (double*)malloc(dimM * sizeof(double));

  for (i = 0; i <= n - 1; i++)
  {
    meth_rp.rp.a[i] = a[i] + c[i] ;
  }


  printf("\n\n we go in the function \n\n");

  info = solve_rp(vec, qqt, &n, &meth_rp, zt, w);


  f7 = fopen("res.dat", "w+");



  for (i = 0; i < n; i++)
  {
    z[i] = zt[i] - c[i] ;
  }


  for (i = 0; i < n; i++)
  {
    fprintf(f7, "%d  %g  %g  %g  %g  \n", i, z[i], w[i], a[i], b[i]);
  }


  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f6);
  fclose(f1);
  fclose(f5);
  fclose(f7);


  free(M);
  free(Mt);
  free(vec);
  free(q);
  free(qqt);
  free(z);
  free(zt);
  free(w);
  free(a);
  free(b);
  free(c);
  free(meth_rp.rp.a);
}

void test2(void)
{
  FILE *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9;
  int i, j, nl, nc, nll, it, info, n = 40, dimM = n, incx, incy;
  double *q, *z, *w, *vec, *a, *b, *c, *qqt, *zt;
  double(*M)[n], (*Mt)[n];
  double qi, Mij, alpha, beta;
  char val[14], vall[14];
  methode meth_rp;
  char trans;

  strcpy(meth_rp.rp.nom_method, "Latin");
  meth_rp.rp.itermax = 50;
  meth_rp.rp.tol = 0.001;
  meth_rp.rp.k_latin = 0.003; //0.00005;



  if ((f1 = fopen("M_relay.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }



  M = malloc(dimM * dimM * sizeof(double));
  Mt = malloc(dimM * dimM * sizeof(double));
  vec = (double*)malloc(dimM * dimM * sizeof(double));


  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);

    /////////////       on met la transpos       ////////////////

    *(*(M + nc - 1) + nl - 1) = Mij;
    *(*(Mt + nl - 1) + nc - 1) = Mij;

    //////////////         fin transpos         ////////////////////

  }

  //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
  for (i = 0; i < dimM; i++)
  {
    for (j = 0; j < dimM; j++)
    {
      vec[j * dimM + i] = M[i][j];

    }
  }
  ////////////////////////////////////////////////////////////////////////




  if ((f2 = fopen("q_relay.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  if ((f5 = fopen("a_relay.dat", "r")) == NULL)
  {
    perror("fopen 5");
    exit(5);
  }


  if ((f6 = fopen("b_relay.dat", "r")) == NULL)
  {
    perror("fopen 6");
    exit(6);
  }



  q = malloc(dimM * sizeof(double));
  qqt = malloc(dimM * sizeof(double));
  z = malloc(dimM * sizeof(double));
  zt = malloc(dimM * sizeof(double));
  w = malloc(dimM * sizeof(double));
  a = malloc(dimM * sizeof(double));
  b = malloc(dimM * sizeof(double));
  c = malloc(dimM * sizeof(double));




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nll);
    fscanf(f2, "%s", vall);
    qi = atof(vall);
    *(q + nll - 1) = qi;
  }


  while (!feof(f5))
  {
    fscanf(f5, "%d", &nll);
    fscanf(f5, "%s", vall);
    qi = atof(vall);
    fprintf(f5, "%d %.14e\n", nll, qi);
    *(a + nll - 1) = qi;
  }


  while (!feof(f6))
  {
    fscanf(f6, "%d", &nll);
    fscanf(f6, "%s", vall);
    qi = atof(vall);
    fprintf(f6, "%d %.14e\n", nll, qi);
    *(b + nll - 1) = qi;
  }


  for (i = 0; i < dimM; i++)
  {
    c[i] = (-b[i] - a[i]) / 2;
    qqt[i] = q[i];
  }


  trans = 'T';
  alpha = 1.;
  beta = 1.;
  incx = 1;
  incy = 1;


  dgemv_(&trans, &n, &n, &alpha, M, &n, c, &incx, &beta, qqt, &incy);

  meth_rp.rp.a = (double*)malloc(dimM * sizeof(double));

  for (i = 0; i <= n - 1; i++)
  {
    meth_rp.rp.a[i] = a[i] + c[i] ;
  }


  printf("\n\n we go in the function \n\n");

  info = solve_rp(vec, qqt, &n, &meth_rp, zt, w);


  f7 = fopen("res.dat", "w+");



  for (i = 0; i < n; i++)
  {
    z[i] = zt[i] - c[i] ;
  }


  for (i = 0; i < n; i++)
  {
    fprintf(f7, "%d  %g  %g  %g  %g  \n", i, z[i], w[i], a[i], b[i]);
  }


  printf("\n\n we go out the function and info is %d\n", info);

  fclose(f2);
  fclose(f6);
  fclose(f1);
  fclose(f5);
  fclose(f7);


  free(M);
  free(Mt);
  free(vec);
  free(q);
  free(qqt);
  free(z);
  free(zt);
  free(w);
  free(a);
  free(b);
  free(c);
  free(meth_rp.rp.a);
}


int main(void)
{
  /*  FILE *f1,*f2,*f3,*f4,*f5,*f6,*f7,*f8,*f9;
    int i,j,nl,nc,nll,it,info,n=40,dimM=n,incx,incy;
    double *q,*z,*w,*vec,*a,*b,*c,*qqt,*zt;
    double (*M)[n],(*Mt)[n];
    double qi,Mij,alpha,beta;
    char val[14],vall[14];
    methode meth_rp;
    char trans;

    //  meth_rp.rp.nom_method="latin";
    strcpy( meth_rp.rp.nom_method, "Gsnl");
    meth_rp.rp.itermax=50;
    meth_rp.rp.tol=0.001;
    meth_rp.rp.k_latin=0.003;//0.00005;



    if ((f1=fopen("M_relay.dat","r"))==NULL){
    perror("fopen 1");
    exit(1);
    }



    M=malloc(dimM*dimM*sizeof(double));
    Mt=malloc(dimM*dimM*sizeof(double));
    vec=(double*)malloc(dimM*dimM*sizeof(double));


    while (!feof(f1)){
      fscanf(f1,"%d",&nl);
      fscanf(f1,"%d",&nc);
      fscanf(f1,"%s",val);
      Mij=atof(val);

   /////////////       on met la transpos       ////////////////

     *(*(M+nc-1)+nl-1)=Mij;
     *(*(Mt+nl-1)+nc-1)=Mij;

   //////////////         fin transpos         ////////////////////

     }

    //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
      for (i=0;i<dimM;i++)
        {for (j=0;j<dimM;j++){
    vec[j*dimM+i]= M[i][j];

        }}
    ////////////////////////////////////////////////////////////////////////




    if ((f2=fopen("q_relay.dat","r"))==NULL){
    perror("fopen 2");
    exit(2);
    }


    if ((f5=fopen("a_relay.dat","r"))==NULL){
    perror("fopen 5");
    exit(5);
    }


    if ((f6=fopen("b_relay.dat","r"))==NULL){
    perror("fopen 6");
    exit(6);
    }



    q=malloc(dimM*sizeof(double));
    qqt=malloc(dimM*sizeof(double));
    z=malloc(dimM*sizeof(double));
    zt=malloc(dimM*sizeof(double));
    w=malloc(dimM*sizeof(double));
    a=malloc(dimM*sizeof(double));
    b=malloc(dimM*sizeof(double));
    c=malloc(dimM*sizeof(double));




    while (!feof(f2)){
      fscanf(f2,"%d",&nll);
      fscanf(f2,"%s",vall);
      qi=atof(vall);
      *(q+nll-1)=qi;
   }


    while (!feof(f5)){
      fscanf(f5,"%d",&nll);
      fscanf(f5,"%s",vall);
      qi=atof(vall);
      fprintf(f5,"%d %.14e\n",nll,qi);
      *(a+nll-1)=qi;
   }


    while (!feof(f6)){
      fscanf(f6,"%d",&nll);
      fscanf(f6,"%s",vall);
      qi=atof(vall);
      fprintf(f6,"%d %.14e\n",nll,qi);
      *(b+nll-1)=qi;
   }


    for(i=0;i<dimM;i++){
      c[i]=(-b[i]-a[i])/2;
      qqt[i]=q[i];
    }


     trans='T';
     alpha=1.;
     beta=1.;
     incx=1;
     incy=1;


    dgemv_(&trans,&n,&n,&alpha,M,&n,c,&incx,&beta,qqt,&incy);

    meth_rp.rp.a=(double*)malloc(dimM*sizeof(double));

    for (i=0;i<=n-1;i++){
      meth_rp.rp.a[i]=a[i]+c[i] ;
    }


    printf("\n\n we go in the function \n\n");

    info=solve_rp (vec,qqt,&n,&meth_rp,zt,w);


    f7=fopen("res.dat","w+");



    for (i=0;i<n;i++){
      z[i]=zt[i]-c[i] ;}


    for(i=0;i<n;i++){
      fprintf(f7,"%d  %g  %g  %g  %g  \n",i,z[i],w[i],a[i],b[i]);
    }


    printf("\n\n we go out the function and info is %d\n",info);

    fclose(f2);fclose(f6);fclose(f1);fclose(f5);fclose(f7);


    free(M);  free(Mt);  free(vec);  free(q);  free(qqt);  free(z);  free(zt);  free(w);  free(a);  free(b);  free(c);  free(meth_rp.rp.a);



  /////////////////////////////////////
  // second test
  ////////////////////////////////////
    strcpy( meth_rp.rp.nom_method, "Latin");
    meth_rp.rp.itermax=50;
    meth_rp.rp.tol=0.001;
    meth_rp.rp.k_latin=0.003;//0.00005;



    if ((f1=fopen("M_relay.dat","r"))==NULL){
    perror("fopen 1");
    exit(1);
    }



    M=malloc(dimM*dimM*sizeof(double));
    Mt=malloc(dimM*dimM*sizeof(double));
    vec=(double*)malloc(dimM*dimM*sizeof(double));


    while (!feof(f1)){
      fscanf(f1,"%d",&nl);
      fscanf(f1,"%d",&nc);
      fscanf(f1,"%s",val);
      Mij=atof(val);

   /////////////       on met la transpos       ////////////////

     *(*(M+nc-1)+nl-1)=Mij;
     *(*(Mt+nl-1)+nc-1)=Mij;

   //////////////         fin transpos         ////////////////////

     }

    //// valeurs du tableau dans vec (compatibilite allocation memoire f90)///
      for (i=0;i<dimM;i++)
        {for (j=0;j<dimM;j++){
    vec[j*dimM+i]= M[i][j];

        }}
    ////////////////////////////////////////////////////////////////////////




    if ((f2=fopen("q_relay.dat","r"))==NULL){
    perror("fopen 2");
    exit(2);
    }


    if ((f5=fopen("a_relay.dat","r"))==NULL){
    perror("fopen 5");
    exit(5);
    }


    if ((f6=fopen("b_relay.dat","r"))==NULL){
    perror("fopen 6");
    exit(6);
    }



    q=malloc(dimM*sizeof(double));
    qqt=malloc(dimM*sizeof(double));
    z=malloc(dimM*sizeof(double));
    zt=malloc(dimM*sizeof(double));
    w=malloc(dimM*sizeof(double));
    a=malloc(dimM*sizeof(double));
    b=malloc(dimM*sizeof(double));
    c=malloc(dimM*sizeof(double));




    while (!feof(f2)){
      fscanf(f2,"%d",&nll);
      fscanf(f2,"%s",vall);
      qi=atof(vall);
      *(q+nll-1)=qi;
   }


    while (!feof(f5)){
      fscanf(f5,"%d",&nll);
      fscanf(f5,"%s",vall);
      qi=atof(vall);
      fprintf(f5,"%d %.14e\n",nll,qi);
      *(a+nll-1)=qi;
   }


    while (!feof(f6)){
      fscanf(f6,"%d",&nll);
      fscanf(f6,"%s",vall);
      qi=atof(vall);
      fprintf(f6,"%d %.14e\n",nll,qi);
      *(b+nll-1)=qi;
   }


    for(i=0;i<dimM;i++){
      c[i]=(-b[i]-a[i])/2;
      qqt[i]=q[i];
    }


     trans='T';
     alpha=1.;
     beta=1.;
     incx=1;
     incy=1;


    dgemv_(&trans,&n,&n,&alpha,M,&n,c,&incx,&beta,qqt,&incy);

    meth_rp.rp.a=(double*)malloc(dimM*sizeof(double));

    for (i=0;i<=n-1;i++){
      meth_rp.rp.a[i]=a[i]+c[i] ;
    }


    printf("\n\n we go in the function \n\n");

    info=solve_rp (vec,qqt,&n,&meth_rp,zt,w);


    f7=fopen("res.dat","w+");



    for (i=0;i<n;i++){
      z[i]=zt[i]-c[i] ;}


    for(i=0;i<n;i++){
      fprintf(f7,"%d  %g  %g  %g  %g  \n",i,z[i],w[i],a[i],b[i]);
    }


    printf("\n\n we go out the function and info is %d\n",info);

    fclose(f2);fclose(f6);fclose(f1);fclose(f5);fclose(f7);


    free(M);  free(Mt);  free(vec);  free(q);  free(qqt);  free(z);  free(zt);  free(w);  free(a);  free(b);  free(c);  free(meth_rp.rp.a);
  */
  test1();
  test2();
  return 1;
}


