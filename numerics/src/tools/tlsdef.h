/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#ifndef TLSDEF_H
#define TLSDEF_H

/*! \file tlsdef.h
 *  \brief definition of thread local variable
 */

#if defined(__GNUC__)
#define DESTRUCTOR_ATTR __attribute__ ((destructor))
#define CONSTRUCTOR_ATTR __attribute__ ((constructor))
#else
#define DESTRUCTOR_ATTR 
#define CONSTRUCTOR_ATTR
#endif


#if !defined(__cplusplus) || !defined(BUILD_AS_CPP)

  #ifdef SWIG
    #define tlsvar

  #elif __STDC_VERSION__ >= 201112L && defined(__STDC_NO_THREADS__) && __STDC_NO_THREADS__ == 0
    #include <threads.h>
    #define tlsvar thread_local
  #else

    #if defined(__GNUC__) || (defined(__ICC) && defined(__linux))
      #define tlsvar __thread
    #elif defined(__ICC) && defined(_WIN32)
      #define tlsvar __declspec(thread)
    #elif defined(SICONOS_ALLOW_GLOBAL)
      #define tlsvar
    #else
      #error "Don't know how to create a thread-local variable"
    #endif
  #endif

#else

  #ifdef SWIG
    #define tlsvar

  #elif __cplusplus >= 201103L
    #define tlsvar thread_local
  #else
    #if defined(__GNUC__) || (defined(__ICC) && defined(__linux))
      #define tlsvar __thread
    #elif defined(_MSC_VER) || (defined(__ICC) && defined(_WIN32))
      #define tlsvar __declspec(thread)
    #elif defined(SICONOS_ALLOW_GLOBAL)
      #define tlsvar
    #else
      #error "Don't know how to create a thread-local variable"
    #endif
  #endif

#endif

#endif
