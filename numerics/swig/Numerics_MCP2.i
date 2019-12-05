%extend MixedComplementarityProblem
{
  CALL_COMPUTE_F(mcp, (void*))

  CALL_COMPUTE_NABLA_F(mcp, (void*))

  MixedComplementarityProblem()
   {
     MixedComplementarityProblem* MCP = mixedComplementarityProblem_new();
     return MCP;
   }

  MixedComplementarityProblem(SN_OBJ_TYPE* n1, SN_OBJ_TYPE* n2)
  {
     MixedComplementarityProblem* MCP = mixedComplementarityProblem_new();

     SWIG_AsVal_int(n1, &MCP->n1);
     SWIG_AsVal_int(n2, &MCP->n2);
     int size =  MCP->n1 +  MCP->n2;

     if (size<1)
     {
       SWIG_Error(SWIG_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       free(MCP);
       return NULL;
     }
     else
     {
     //TODO implement different types of matrices
       MCP->nabla_Fmcp = NM_create(NM_DENSE, size, size);
     }
     return MCP;
  }


#ifdef SWIGPYTHON
  MixedComplementarityProblem(SN_OBJ_TYPE* n1, SN_OBJ_TYPE* n2, SN_OBJ_TYPE* py_compute)
  {
     MixedComplementarityProblem* MCP = mixedComplementarityProblem_new();

     MCP->compute_Fmcp = &call_py_compute_Fmcp;
     MCP->compute_nabla_Fmcp = &call_py_compute_nabla_Fmcp;
     SWIG_AsVal_int(n1, &MCP->n1);
     SWIG_AsVal_int(n2, &MCP->n2);
     int size =  MCP->n1 +  MCP->n2;

     if (size<1)
     {
       SWIG_Error(SWIG_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       free(MCP);
       return NULL;
     }
     else
     {
       MCP->nabla_Fmcp = NM_create(NM_DENSE, size, size);
     }

     SN_OBJ_TYPE* method_compute_F = NULL;
     if (PyObject_HasAttrString(py_compute, "compute_F")) method_compute_F = PyObject_GetAttrString(py_compute, "compute_F");
     SN_OBJ_TYPE* method_compute_nabla_F = NULL;
     if (PyObject_HasAttrString(py_compute, "compute_nabla_F")) method_compute_nabla_F = PyObject_GetAttrString(py_compute, "compute_nabla_F");
     if (PyCallable_Check(method_compute_F) && PyCallable_Check(method_compute_nabla_F))
     {
       MCP->env = (void*) malloc(sizeof(class_env_python));
       class_env_python* mcp_env_python = (class_env_python*) MCP->env;
       mcp_env_python->id = ENV_IS_PYTHON_CLASS;
       mcp_env_python->class_object = py_compute;
       target_mem_mgmt_instr(method_compute_F);
       target_mem_mgmt_instr(method_compute_nabla_F);
     }
     else
     {
       target_mem_mgmtX_instr(method_compute_F);
       target_mem_mgmtX_instr(method_compute_nabla_F);
       SWIG_Error(SWIG_TypeError, "argument 2 must be have a method compute_F and a method compute_nabla_F");
       NM_clear(MCP->nabla_Fmcp);
       free(MCP->nabla_Fmcp);
       free(MCP);
       return NULL;
     }

     return MCP;
   }
#endif /* SWIGPYTHON */

  MixedComplementarityProblem(SN_OBJ_TYPE* n1, SN_OBJ_TYPE* n2, SN_OBJ_TYPE* compute_F, SN_OBJ_TYPE* compute_nabla_F)
  {
     MixedComplementarityProblem* MCP = mixedComplementarityProblem_new();

     SWIG_AsVal_int(n1, &MCP->n1);
     SWIG_AsVal_int(n2, &MCP->n2);
     int size =  MCP->n1 +  MCP->n2;

     if (size<1)
     {
       SWIG_Error(SWIG_RuntimeError, "sizeEqualities + sizeInequalities has to be positive");
       free(MCP);
       return NULL;
     }
     else
     {
       MCP->nabla_Fmcp = NM_create(NM_DENSE, size, size);
     }

     check_save_target_fn(compute_F, MCP->env, env_compute_function, MixedComplementarityProblem_call_compute_F, MCP->compute_Fmcp, 2);
     check_save_target_fn(compute_nabla_F, MCP->env, env_compute_jacobian, MixedComplementarityProblem_call_compute_nabla_F, MCP->compute_nabla_Fmcp, 3);

     return MCP;
   }

#ifdef SWIGPYTHON
  void set_compute_F_and_nabla_F_as_C_functions(SN_OBJ_TYPE* lib_name, SN_OBJ_TYPE* compute_F_name, SN_OBJ_TYPE* compute_nabla_F_name)
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
      SWIG_Error(SWIG_TypeError, "All arguments should be strings");
    }
  }
#endif /* SWIGPYTHON */

    SN_OBJ_TYPE* get_env_as_long(void)
    {
      return SWIG_From_long((uintptr_t)&$self->env);
    }

  ~MixedComplementarityProblem()
  {
    if ($self->env)
    {
      if(((env_target_lang*)$self->env)->id > 0)
      {
        free($self->env);
        $self->env = NULL;
      }
    }
    mixedComplementarityProblem_free($self);
  }
};


