
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

static void call_py_compute_nabla_Fvi(void *env, int n, double* z, NumericsMatrix* nabla_F);
static void call_py_compute_Fvi(void *env, int n, double* z, double* F);

static void call_py_compute_nabla_Fncp(void *env, int n, double* z, NumericsMatrix* nabla_F);
static void call_py_compute_Fncp(void *env, int n, double* z, double* F);

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

//////////////////////////////////////////////////////////////////////////////
// START copy from npy_3kcompat.h
// We copy those things since
//   - there is no guarantee that those functions are going to remain around
//   - MSVC chokes if we include npy_3kcompat.h and I don't want to investigate
// 
// --xhub
//////////////////////////////////////////////////////////////////////////////

#if PY_MAJOR_VERSION >= 3

#define PyString_Type PyBytes_Type
#define PyStringObject PyBytesObject
#define PyString_FromStringAndSize PyBytes_FromStringAndSize
#define PyString_AsStringAndSize PyBytes_AsStringAndSize
#define PyString_FromFormat PyBytes_FromFormat
#define PyString_Concat PyBytes_Concat
#define PyString_ConcatAndDel PyBytes_ConcatAndDel
#define PyString_GET_SIZE PyBytes_GET_SIZE

#define PyUString_Type PyUnicode_Type
#define PyUString_Check PyUnicode_Check
#define PyUStringObject PyUnicodeObject
#define PyUString_FromString PyUnicode_FromString
#define PyUString_FromStringAndSize PyUnicode_FromStringAndSize
#define PyUString_FromFormat PyUnicode_FromFormat
#define PyUString_Concat PyUnicode_Concat2
#define PyUString_ConcatAndDel PyUnicode_ConcatAndDel
#define PyUString_GET_SIZE PyUnicode_GET_SIZE
#define PyUString_Size PyUnicode_Size
#define PyUString_InternFromString PyUnicode_InternFromString
#define PyUString_Format PyUnicode_Format

#else

#define PyBytes_Type PyString_Type
#define PyBytes_Check PyString_Check
#define PyBytesObject PyStringObject
#define PyBytes_FromString PyString_FromString
#define PyBytes_FromStringAndSize PyString_FromStringAndSize
#define PyBytes_AS_STRING PyString_AS_STRING
#define PyBytes_AsStringAndSize PyString_AsStringAndSize
#define PyBytes_FromFormat PyString_FromFormat
#define PyBytes_Concat PyString_Concat
#define PyBytes_ConcatAndDel PyString_ConcatAndDel
#define PyBytes_AsString PyString_AsString
#define PyBytes_GET_SIZE PyString_GET_SIZE
#define PyBytes_Size PyString_Size

#define PyUString_Type PyString_Type
#define PyUString_Check PyString_Check
#define PyUStringObject PyStringObject
#define PyUString_FromString PyString_FromString
#define PyUString_FromStringAndSize PyString_FromStringAndSize
#define PyUString_FromFormat PyString_FromFormat
#define PyUString_Concat PyString_Concat
#define PyUString_ConcatAndDel PyString_ConcatAndDel
#define PyUString_GET_SIZE PyString_GET_SIZE
#define PyUString_Size PyString_Size
#define PyUString_InternFromString PyString_InternFromString
#define PyUString_Format PyString_Format

#endif /* NPY_PY3K */


static inline void
PyUnicode_ConcatAndDel(PyObject **left, PyObject *right)
{
    PyObject *newobj;
    newobj = PyUnicode_Concat(*left, right);
    Py_DECREF(*left);
    Py_DECREF(right);
    *left = newobj;
}

//////////////////////////////////////////////////////////////////////////////
// END copy from npy_3kcompat.h
//////////////////////////////////////////////////////////////////////////////

