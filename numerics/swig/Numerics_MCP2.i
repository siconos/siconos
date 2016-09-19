%extend MixedComplementarityProblem2
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

  MixedComplementarityProblem2(PyObject* n1, PyObject* n2)
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
     //TODO implement different types of matrices
       MCP->nabla_Fmcp = createNumericsMatrix(NM_DENSE, size, size);
     }
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
       MCP->nabla_Fmcp = createNumericsMatrix(NM_DENSE, size, size);
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
       freeNumericsMatrix(MCP->nabla_Fmcp);
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
       MCP->nabla_Fmcp = createNumericsMatrix(NM_DENSE, size, size);
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
       freeNumericsMatrix(MCP->nabla_Fmcp);
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
       freeNumericsMatrix(MCP->nabla_Fmcp);
       free(MCP->nabla_Fmcp);
       free(MCP->env);
       free(MCP);
       PyErr_PrintEx(0);
       exit(1);
     }

     return MCP;
   }

  void set_compute_F_and_nabla_F_as_C_functions(PyObject* lib_name, PyObject* compute_F_name, PyObject* compute_nabla_F_name)
  {
%#if PY_MAJOR_VERSION < 3
    if(PyString_Check(lib_name) && PyString_Check(compute_F_name) && PyString_Check(compute_nabla_F_name))
%#else
    if(PyUnicode_Check(lib_name) && PyUnicode_Check(compute_F_name) && PyUnicode_Check(compute_nabla_F_name))
%#endif
    {
    void* p_compute_F;
    void* p_compute_nabla_F;

    // TODO: save this lib_handle somewhere and close it !
    get_c_functions(lib_name, compute_F_name, compute_nabla_F_name, &p_compute_F, &p_compute_nabla_F);

    $self->compute_Fmcp = (ptrFunctionMCP2)p_compute_F;
    $self->compute_nabla_Fmcp = (ptrFunctionMCP_nabla)p_compute_nabla_F;

    }
    else
    {
      PyErr_SetString(PyExc_TypeError, "All arguments should be strings");
    }
  }

    PyObject* get_env_as_long(void)
    {
      return PyInt_FromLong((uintptr_t)&$self->env);
    }

  ~MixedComplementarityProblem2()
  {
    if ($self->nabla_Fmcp)
    {
      freeNumericsMatrix($self->nabla_Fmcp);
      free($self->nabla_Fmcp);
    }
    if ($self->env)
    {
      if(((env_python*)$self->env)->id > 0)
      {
        free($self->env);
      }
    }
    free($self);
  }
};


