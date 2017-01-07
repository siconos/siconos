/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*!\file sn_error_handling.h
 * \brief error handling functions and data structures*/

#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#include "SiconosConfig.h"
#include "tlsdef.h"
#include <setjmp.h>
#include <stdbool.h>

extern tlsvar jmp_buf external_jmp_buf;
extern tlsvar jmp_buf internal_jmp_buf;
extern tlsvar bool external_jmp_buf_used;
extern tlsvar bool internal_jmp_buf_used;

#define SN_SETJMP_INTERNAL_START setjmp(internal_jmp_buf), internal_jmp_buf_used = true
#define SN_STEJMP_INTERNAL_STOP internal_jmp_buf_used = false

#define SN_SETJMP_EXTERNAL_START setjmp(*sn_get_jmp_buf())
#define SN_SETJMP_EXTERNAL_STOP sn_release_jmp_buf();

typedef enum { SN_NO_ERROR, SN_MEMORY_ALLOC_ERROR, SN_UNSUPPORTED_LINALG_OP, SN_PROBLEM_NOT_PROCESSABLE, SN_UNKOWN_ERROR, SN_NOT_COMPILED_ERROR } SN_ERROR_T;

#if defined(__cplusplus) && !defined (BUILD_AS_CPP)
extern "C"
{
#endif


  jmp_buf* sn_get_jmp_buf(void);
  void sn_release_jmp_buf(void);
  void sn_fatal_error(SN_ERROR_T code, const char* msg);
  const char* sn_fatal_error_msg(void);

#if defined(__cplusplus) && !defined (BUILD_AS_CPP)
}
#endif

#endif /* ERROR_HANDLING_H  */
