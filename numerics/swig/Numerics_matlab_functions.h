#include "Numerics_common_structs.h"

typedef struct {
  int id;
  void* env_compute_function;
  void* env_compute_jacobian;
  void* env_compute_function_str;
  void* env_compute_jacobian_str;
} functions_env_matlab;
