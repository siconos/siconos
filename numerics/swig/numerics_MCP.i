// MCP
%{
#include "MixedComplementarityProblem.h"
#include "MCP_Solvers.h"
#include "MCP_cst.h"
  %}

%include "MixedComplementarityProblem.h"
%include "MCP_Solvers.h"
%include "MCP_cst.h"

%extend MixedComplementarityProblem
{



  MixedComplementarityProblem()
   {
     MixedComplementarityProblem* MCP;
     MCP =  (MixedComplementarityProblem *) malloc(sizeof(MixedComplementarityProblem));
     MCP->Fmcp=NULL;
     MCP->nablaFmcp=NULL;
     MCP->computeFmcp=NULL;
     MCP->computeNablaFmcp=NULL;
     return MCP;
   }

  void set_computeFmcp(PyObject *o)
  {
    set_my_callback_Fmcp(o);
    $self->computeFmcp = (my_call_to_callback_Fmcp);
  }

  void set_computeNablaFmcp(PyObject *o)
  {

    set_my_callback_NablaFmcp(o);
    $self->computeNablaFmcp = (my_call_to_callback_NablaFmcp);
  }

  void test_call_to_callback()
  {
    printf("I am in test_call_to_callback()\n");

    int size =   $self->sizeEqualities +  $self->sizeInequalities;

    double * z = (double *)malloc(size*sizeof(double));
    double * F = (double *)malloc(size*sizeof(double));
    double * nablaF = (double *)malloc(size*size*sizeof(double));

    for (int i=0; i < size; i++) z[i]=i;
    printf("Input \n");
    for (int i=0; i < size; i++) printf("z[%i] = %lf\t", i, z[i]);
    printf("\n");
    $self->computeFmcp(size,z,F);
    if  (!PyErr_Occurred())
    {
      $self->computeNablaFmcp(size,z,nablaF);
    }
    printf("Output \n");
    for (int i=0; i < size; i++) printf("F[%i] =  %lf\t", i, F[i]);
    printf("\n");
    for (int i=0; i < size*size; i++) printf("nablaF[%i] =  %lf\t", i, nablaF[i]);

    printf("\n");
    free(z);
    free(F);
    free(nablaF);



    printf("I leave test_call_to_callback()\n");
  }



  MixedComplementarityProblem(PyObject *sizeEq, PyObject *sizeIneq, PyObject *o1, PyObject *o2)
  {
     MixedComplementarityProblem* MCP;
     MCP =  (MixedComplementarityProblem *) malloc(sizeof(MixedComplementarityProblem));

     MCP->sizeEqualities = (int) PyInt_AsLong(sizeEq);
     MCP->sizeInequalities = (int) PyInt_AsLong(sizeIneq);
     int size =  MCP->sizeEqualities +  MCP->sizeInequalities;

     if (size<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       MCP->Fmcp = NULL;
       MCP->nablaFmcp = NULL;
       freeMixedComplementarityProblem(MCP);
       return NULL;
     }
     else
     {
       MCP->Fmcp = (double *) malloc(size*sizeof(double));
       MCP->nablaFmcp = (double *) malloc(size*size*sizeof(double));
     }

     if (PyCallable_Check(o1))
     {
       set_my_callback_Fmcp(o1);
       MCP->computeFmcp = (my_call_to_callback_Fmcp);
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 3 must be callable");
       free(MCP->Fmcp);
       free(MCP->nablaFmcp);
       MCP->Fmcp = NULL;
       MCP->nablaFmcp = NULL;
       freeMixedComplementarityProblem(MCP);
       return NULL;
     }


     if (PyCallable_Check(o2))
     {
       set_my_callback_NablaFmcp(o2);
       MCP->computeNablaFmcp = (my_call_to_callback_NablaFmcp);
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 4 must be callable");
       free(MCP->Fmcp);
       free(MCP->nablaFmcp);
       MCP->Fmcp = NULL;
       MCP->nablaFmcp = NULL;
       freeMixedComplementarityProblem(MCP);
       return NULL;
     }

     return MCP;
   }

  ~MixedComplementarityProblem()
  {
    free($self->Fmcp);
    free($self->nablaFmcp);
    $self->Fmcp = NULL;
    $self->nablaFmcp = NULL;
    freeMixedComplementarityProblem($self);
  }
};



