
%extend MixedComplementarityProblem2_
{
  MixedComplementarityProblem2()
   {
     MixedComplementarityProblem2* MCP;
     MCP = (MixedComplementarityProblem2 *) malloc(sizeof(MixedComplementarityProblem2));
     MCP->nabla_Fmcp = NULL;
     MCP->compute_Fmcp = &call_py_compute_Fmcp;
     MCP->compute_nabla_Fmcp = &call_py_compute_nabla_Fmcp;
     MCP->env = NULL;

     return MCP;
   }

  MixedComplementarityProblem2(PyObject* n1, PyObject* n2, PyObject* py_compute_class)
  {
     MixedComplementarityProblem2* MCP;
     MCP =  (MixedComplementarityProblem2 *) malloc(sizeof(MixedComplementarityProblem2));

     MCP->compute_Fmcp = &call_py_compute_Fmcp;
     MCP->compute_nabla_Fmcp = &call_py_compute_nabla_Fmcp;
     MCP->n1 = (int) PyInt_AsLong(n1);
     MCP->n2 = (int) PyInt_AsLong(n2);
     int size =  MCP->n1 +  MCP->n2;

     if (size<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       free(MCP);
       return NULL;
     }
     else
     {
       MCP->nabla_Fmcp = (double *) calloc(size*size, sizeof(double));
     }

     PyObject* method_compute_Fmcp = PyObject_GetAttrString(py_compute_class, "compute_Fmcp");
     PyObject* method_compute_nabla_Fmcp = PyObject_GetAttrString(py_compute_class, "compute_nabla_Fmcp");
     if (PyCallable_Check(method_compute_Fmcp) && PyCallable_Check(method_compute_nabla_Fmcp))
     {
       MCP->env = (void*) malloc(sizeof(class_env_python));
       class_env_python* mcp_env_python = (class_env_python*) MCP->env;
       mcp_env_python->id = ENV_IS_PYTHON_CLASS;
       mcp_env_python->class_object = py_compute_class;
       Py_DECREF(method_compute_Fmcp);
       Py_DECREF(method_compute_nabla_Fmcp);
     }
     else
     {
       Py_XDECREF(method_compute_Fmcp);
       Py_XDECREF(method_compute_nabla_Fmcp);
       PyErr_SetString(PyExc_TypeError, "argument 2 must be have a method compute_Fmcp and a method compute_nabla_Fmcp");
       free(MCP->nabla_Fmcp);
       free(MCP);
       PyErr_PrintEx(0);
       exit(1);
     }

     return MCP;
   }

  MixedComplementarityProblem2(PyObject* n1, PyObject* n2, PyObject* py_compute_Fmcp, PyObject* py_compute_nabla_Fmcp)
  {
     MixedComplementarityProblem2* MCP;
     MCP =  (MixedComplementarityProblem2 *) malloc(sizeof(MixedComplementarityProblem2));

     MCP->compute_Fmcp = &call_py_compute_Fmcp;
     MCP->compute_nabla_Fmcp = &call_py_compute_nabla_Fmcp;
     MCP->n1 = (int) PyInt_AsLong(n1);
     MCP->n2 = (int) PyInt_AsLong(n2);
     int size =  MCP->n1 +  MCP->n2;

     if (size<1)
     {
       PyErr_SetString(PyExc_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       free(MCP);
       return NULL;
     }
     else
     {
       MCP->nabla_Fmcp = (double *) malloc(size*size*sizeof(double));
     }

     MCP->env = (void*) malloc(sizeof(functions_env_python));
     functions_env_python* mcp_env_python = (functions_env_python*) MCP->env;
     mcp_env_python->id = ENV_IS_PYTHON_FUNCTIONS;

     if (PyCallable_Check(py_compute_Fmcp)) 
     {
       mcp_env_python->env_compute_function = py_compute_Fmcp;
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 3 must be callable");
       free(MCP->nabla_Fmcp);
       free(MCP->env);
       free(MCP);
       PyErr_PrintEx(0);
       exit(1);
     }


     if (PyCallable_Check(py_compute_nabla_Fmcp))
     {
       mcp_env_python->env_compute_jacobian = py_compute_nabla_Fmcp;
     }
     else
     {
       PyErr_SetString(PyExc_TypeError, "argument 4 must be callable");
       free(MCP->nabla_Fmcp);
       free(MCP->env);
       free(MCP);
       PyErr_PrintEx(0);
       exit(1);
     }

     return MCP;
   }
  ~MixedComplementarityProblem2()
  {
    free($self->nabla_Fmcp);
    free($self->env);
    free($self);
  }
};


