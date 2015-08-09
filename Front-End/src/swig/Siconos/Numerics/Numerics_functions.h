
// Heritage in C is wonderful

typedef struct {
  int id;
} env_python;

typedef struct {
  int id;
  PyObject* class_object;
} class_env_python;

typedef struct {
  int id;
  PyObject* env_compute_function;
  PyObject* env_compute_jacobian;
} functions_env_python;

typedef struct {
  int id;
  PyObject* env_compute_function;
  PyObject* env_compute_projection_on_set;
} functions_env_python_with_proj;


#define ENV_IS_PYTHON_CLASS 1
#define ENV_IS_PYTHON_FUNCTIONS 2
#define ENV_IS_PYTHON_FUNCTIONS_WITH_PROJ 3
#define ENV_IS_C_STRUCT -1

static void call_py_compute_nabla_F(void *env, int n, double* z, NumericsMatrix* nabla_F);
static void call_py_compute_F(void *env, int n, double* z, double* F);

static void call_py_compute_nabla_Fmcp(void *env, int n1, int n2, double* z, NumericsMatrix* nabla_Fmcp);
static void call_py_compute_Fmcp(void *env, int n1, int n2, double* z, double* Fmcp);

#define PYTHON_NAME_CALL_TYPE_2 PyString_FromString
#define PYTHON_NAME_CALL_TYPE_3 PyUnicode_FromString
#define PYTHON_NAME_CALL_TYPE_4 I_see_no_future

#define CONC(A,B) A##_##B

#define GET_PYTHON_NAME_CALL(X) CONC(PYTHON_NAME_CALL_TYPE, X)

#define PY_CALL_METHOD_OR_FUNCTION(ENV_STRUCT, METHOD_NAME, FIELD_FUNCTION, ...) \
PyObject* py_out = NULL;\
switch (((env_python*) ENV_STRUCT)->id)\
{\
  case ENV_IS_PYTHON_CLASS:\
  {\
    py_out = PyObject_CallMethodObjArgs(((class_env_python*) ENV_STRUCT)->class_object, GET_PYTHON_NAME_CALL(PY_MAJOR_VERSION)((char *)METHOD_NAME), __VA_ARGS__, NULL);\
    break;\
  }\
  case ENV_IS_PYTHON_FUNCTIONS:\
  {\
    py_out = PyObject_CallFunctionObjArgs(((functions_env_python*) ENV_STRUCT)->FIELD_FUNCTION, __VA_ARGS__, NULL);\
    break;\
  }\
  default:\
  {\
    PyErr_SetString(PyExc_TypeError, "Unknown environment type");\
  }\
}\
if (py_out == NULL)\
{\
  PyErr_PrintEx(0);\
  exit(1);\
}\
Py_XDECREF(py_out);

