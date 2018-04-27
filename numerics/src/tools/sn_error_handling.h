/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#define SN_SETJMP_INTERNAL_START setjmp(*sn_get_internal_jmp_buf())
#define SN_SETJMP_INTERNAL_STOP sn_release_internal_jmp_buf();

#define SN_SETJMP_EXTERNAL_START setjmp(*sn_get_jmp_buf())
#define SN_SETJMP_EXTERNAL_STOP sn_release_jmp_buf();

typedef enum { SN_NO_ERROR, SN_MEMORY_ALLOC_ERROR, SN_UNSUPPORTED_LINALG_OP, SN_PROBLEM_NOT_PROCESSABLE, SN_UNKOWN_ERROR, SN_NOT_COMPILED_ERROR } SN_ERROR_T;

#if defined(__cplusplus) && !defined (BUILD_AS_CPP)
extern "C"
{
#endif


  /* Get the external jmp buffer and mark it as used
   * \return the external jmp buffer
   */
  jmp_buf* sn_get_jmp_buf(void);

  /** Release the internal jmp buffer: this indicates that it is no longer in
   * use and that there should be no longjmp() call to this part of the stack.
   * The user should call this function whenever the call to a numerics
   * function has been successful*/
  void sn_release_jmp_buf(void);

  /* Get the internal jmp buffer and mark it as used
   * \warning this function is ment to be called inside the numerics library.
   * To use the exception handler from an external library/executable, use
   * sn_get_jmp_buf()
   * \return the internal jmp buffer
   */
  jmp_buf* sn_get_internal_jmp_buf(void);

  /** Release the internal jmp buffer: this indicates that it is no longer in
   * use and that there should be no longjmp() call to this part of the stack.
   * The user should call this function whenever the call to a numerics
   * function has been successful.
   * \warning This should not be called outside the numerics library. Use
   * sn_release_jmp_buf() instead when calling from another library/executable
   */
  void sn_release_internal_jmp_buf(void);

  /* Function to call whenever a fatal error occured. This function either call
   * longjmp if setjmp has been called previously and is still active. If not,
   * it calls abort().
   * \param code error code
   * \param msn error message
   */
  void sn_fatal_error(SN_ERROR_T code, const char* msg);

  /* Get the last error message
   * \return the error message
   */
  const char* sn_fatal_error_msg(void);

#if defined(__cplusplus) && !defined (BUILD_AS_CPP)
}
#endif

#endif /* ERROR_HANDLING_H  */
