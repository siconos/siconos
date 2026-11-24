// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "sn_error_handling.h"
#ifndef __cplusplus
#include <stdbool.h>  // for false, bool, true
#endif
#include <stdlib.h>  // for NULL, abort, size_t

#include "numerics_errors.h"
#include "tlsdef.h"  // for tlsvar

tlsvar jmp_buf internal_jmp_buf;
tlsvar jmp_buf external_jmp_buf;

tlsvar bool external_jmp_buf_used = false;
tlsvar bool internal_jmp_buf_used = false;

tlsvar const char* internal_jmp_buf_err = NULL;
tlsvar const char* external_jmp_buf_err = NULL;

typedef void (*external_fault_handler_t)(size_t, const char*);
tlsvar external_fault_handler_t external_fault_handler = NULL;

jmp_buf* sn_get_jmp_buf(void) {
  external_jmp_buf_used = true;
  return &external_jmp_buf;
}

void sn_release_jmp_buf(void) { external_jmp_buf_used = false; }

jmp_buf* sn_get_internal_jmp_buf(void) {
  internal_jmp_buf_used = true;
  return &internal_jmp_buf;
}

void sn_release_internal_jmp_buf(void) { internal_jmp_buf_used = false; }

void sn_fatal_error(SN_ERROR_T code, const char* msg) {
  if (external_fault_handler) {
    (*external_fault_handler)(code, msg);
  }

  if (internal_jmp_buf_used) {
    internal_jmp_buf_used = false;
    internal_jmp_buf_err = msg;
    longjmp(internal_jmp_buf, code);
  } else if (external_jmp_buf_used) {
    external_jmp_buf_used = false;
    external_jmp_buf_err = msg;
    longjmp(external_jmp_buf, code);
  } else {
    abort();
  }
}

const char* sn_fatal_error_msg(void) {
  const char* err_msg = NULL;
  if (internal_jmp_buf_err) {
    err_msg = internal_jmp_buf_err;
    internal_jmp_buf_err = NULL;
  } else if (external_jmp_buf_err) {
    err_msg = external_jmp_buf_err;
    external_jmp_buf_err = NULL;
  }
  return err_msg;
}
