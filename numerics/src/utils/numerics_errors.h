/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*!\file numerics_errors.h
 * \brief Standardized error codes and error handling macros
 *
 * This header provides consistent error handling across all numerics solvers,
 * replacing the current mix of return codes, assertions, and error macros.
 */

#ifndef NUMERICS_ERRORS_H
#define NUMERICS_ERRORS_H

#include <stdio.h> /* for fprintf */

#ifdef __cplusplus
extern "C" {
#endif

/* ===========================================================================
 * Error Codes
 * =========================================================================== */

typedef enum {
    /* Success */
    NUMERICS_OK = 0,
    
    /* General errors (1-19) */
    NUMERICS_ERR_NULL_POINTER = -1,
    NUMERICS_ERR_INVALID_ARGUMENT = -2,
    NUMERICS_ERR_INVALID_OPTION = -3,
    NUMERICS_ERR_INVALID_SOLVER = -4,
    NUMERICS_ERR_NOT_IMPLEMENTED = -5,
    NUMERICS_ERR_MEMORY = -6,
    NUMERICS_ERR_FILE_IO = -7,
    
    /* Solver-specific errors (20-39) */
    NUMERICS_ERR_MAX_ITER = -20,
    NUMERICS_ERR_DIVERGENCE = -21,
    NUMERICS_ERR_STAGNATION = -22,
    NUMERICS_ERR_LINEAR_SOLVER = -23,
    NUMERICS_ERR_LOCAL_SOLVER = -24,
    NUMERICS_ERR_PROJECTION = -25,
    
    /* Problem-specific errors (40-59) */
    NUMERICS_ERR_INFEASIBLE = -40,
    NUMERICS_ERR_UNBOUNDED = -41,
    NUMERICS_ERR_ILL_CONDITIONED = -42,
    NUMERICS_ERR_SINGULAR_MATRIX = -43,
    
    /* Unknown error */
    NUMERICS_ERR_UNKNOWN = -99
} NumericsError;

/* ===========================================================================
 * Error Messages
 * =========================================================================== */

static inline const char* numerics_error_string(NumericsError err) {
    switch (err) {
        case NUMERICS_OK: return "Success";
        case NUMERICS_ERR_NULL_POINTER: return "Null pointer";
        case NUMERICS_ERR_INVALID_ARGUMENT: return "Invalid argument";
        case NUMERICS_ERR_INVALID_OPTION: return "Invalid option";
        case NUMERICS_ERR_INVALID_SOLVER: return "Invalid solver";
        case NUMERICS_ERR_NOT_IMPLEMENTED: return "Not implemented";
        case NUMERICS_ERR_MEMORY: return "Memory allocation failed";
        case NUMERICS_ERR_FILE_IO: return "File I/O error";
        case NUMERICS_ERR_MAX_ITER: return "Maximum iterations reached";
        case NUMERICS_ERR_DIVERGENCE: return "Solver diverged";
        case NUMERICS_ERR_STAGNATION: return "Solver stagnated";
        case NUMERICS_ERR_LINEAR_SOLVER: return "Linear solver failed";
        case NUMERICS_ERR_LOCAL_SOLVER: return "Local solver failed";
        case NUMERICS_ERR_PROJECTION: return "Projection failed";
        case NUMERICS_ERR_INFEASIBLE: return "Problem infeasible";
        case NUMERICS_ERR_UNBOUNDED: return "Problem unbounded";
        case NUMERICS_ERR_ILL_CONDITIONED: return "Ill-conditioned problem";
        case NUMERICS_ERR_SINGULAR_MATRIX: return "Singular matrix";
        default: return "Unknown error";
    }
}

/* ===========================================================================
 * Error Checking Macros
 * =========================================================================== */

/** Check condition, return error if false */
#define NUMERICS_CHECK(cond, err, ...) \
    do { \
        if (!(cond)) { \
            numerics_log_error(__FILE__, __LINE__, __func__, \
                               numerics_error_string(err), __VA_ARGS__); \
            return err; \
        } \
    } while(0)

/** Check pointer is not NULL */
#define NUMERICS_CHECK_NULL(ptr) \
    NUMERICS_CHECK((ptr) != NULL, NUMERICS_ERR_NULL_POINTER, \
                   #ptr " is NULL")

/** Check solver ID is valid (solver exists in registry) */
#define NUMERICS_CHECK_SOLVER(id) \
    NUMERICS_CHECK(solver_registry_exists(id), NUMERICS_ERR_INVALID_SOLVER, \
                   "Solver ID %d not registered", (int)(id))

/** Check option index is within bounds */
#define NUMERICS_CHECK_IPARAM(options, idx) \
    NUMERICS_CHECK((idx) >= 0 && (idx) < (options)->iSize, \
                   NUMERICS_ERR_INVALID_OPTION, \
                   "iparam index %d out of bounds (size=%d)", \
                   (int)(idx), (options)->iSize)

#define NUMERICS_CHECK_DPARAM(options, idx) \
    NUMERICS_CHECK((idx) >= 0 && (idx) < (options)->dSize, \
                   NUMERICS_ERR_INVALID_OPTION, \
                   "dparam index %d out of bounds (size=%d)", \
                   (int)(idx), (options)->dSize)

/** Check solver converged (error < tolerance) */
#define NUMERICS_CHECK_CONVERGENCE(error, tol) \
    do { \
        if ((error) > (tol)) { \
            numerics_log_error(__FILE__, __LINE__, __func__, \
                               numerics_error_string(NUMERICS_ERR_DIVERGENCE), \
                               "Error %e > tolerance %e", (error), (tol)); \
            return NUMERICS_ERR_DIVERGENCE; \
        } \
    } while(0)

/** Check iteration limit */
#define NUMERICS_CHECK_ITER(iter, max_iter) \
    do { \
        if ((iter) >= (max_iter)) { \
            numerics_log_error(__FILE__, __LINE__, __func__, \
                               numerics_error_string(NUMERICS_ERR_MAX_ITER), \
                               "Iteration %d >= max %d", (int)(iter), (int)(max_iter)); \
            return NUMERICS_ERR_MAX_ITER; \
        } \
    } while(0)

/* ===========================================================================
 * Logging Functions
 * =========================================================================== */

void numerics_log_error(const char* file, int line, const char* func,
                        const char* error_type, const char* fmt, ...);

void numerics_log_warning(const char* file, int line, const char* func,
                          const char* fmt, ...);

void numerics_log_info(const char* file, int line, const char* func,
                       const char* fmt, ...);

/* ===========================================================================
 * Standardized Error Handling Macros
 * =========================================================================== */

/** Standardized NULL pointer check - returns NUMERICS_ERR_NULL_POINTER */
#define CHECK_NULL(ptr) \
    do { \
        if ((ptr) == NULL) { \
            fprintf(stderr, "[ERROR] %s:%d (%s) - Null pointer: " #ptr "\n", \
                    __FILE__, __LINE__, __func__); \
            return NUMERICS_ERR_NULL_POINTER; \
        } \
    } while(0)


/** Standardized NULL pointer check - returns NUMERICS_ERR_NULL_POINTER */
#define CHECK_ZERO(value) \
    do { \
        if ((value) == 0) { \
            fprintf(stderr, "[ERROR] %s:%d (%s) - zero value: " #value "\n", \
                    __FILE__, __LINE__, __func__); \
            return NUMERICS_ERR_NULL_POINTER; \
        } \
    } while(0)

/** Check condition and return specific error code */
#define CHECK_COND(cond, err_code, msg) \
    do { \
        if (!(cond)) { \
            fprintf(stderr, "[ERROR] %s:%d (%s) - %s\n", \
                    __FILE__, __LINE__, __func__, msg); \
            return err_code; \
        } \
    } while(0)

/** Check condition and return NUMERICS_ERR_INVALID_ARGUMENT */
#define CHECK_ARG(cond, msg) CHECK_COND(cond, NUMERICS_ERR_INVALID_ARGUMENT, msg)

/** Check solver options validity */
#define CHECK_OPTIONS(opts) \
    do { \
        CHECK_NULL(opts); \
        CHECK_COND((opts)->iparam != NULL, NUMERICS_ERR_INVALID_OPTION, "iparam is NULL"); \
        CHECK_COND((opts)->dparam != NULL, NUMERICS_ERR_INVALID_OPTION, "dparam is NULL"); \
    } while(0)

/** Check problem dimensions are valid */
#define CHECK_DIMENSION(dim, expected) \
    CHECK_COND((dim) == (expected), NUMERICS_ERR_INVALID_ARGUMENT, \
               "Invalid dimension: expected " #expected ", got " #dim)

/** Check matrix is not NULL and has valid dimensions */
#define CHECK_MATRIX(M) \
    do { \
        CHECK_NULL(M); \
        CHECK_COND((M)->size0 > 0 && (M)->size1 > 0, NUMERICS_ERR_INVALID_ARGUMENT, \
                   "Matrix has invalid dimensions"); \
    } while(0)

/** Standardized iteration limit check with error logging */
#define CHECK_MAX_ITER(iter, max_iter) \
    do { \
        if ((iter) >= (max_iter)) { \
            fprintf(stderr, "[ERROR] %s:%d (%s) - Maximum iterations reached: %d >= %d\n", \
                    __FILE__, __LINE__, __func__, (int)(iter), (int)(max_iter)); \
            return NUMERICS_ERR_MAX_ITER; \
        } \
    } while(0)

/** Check memory allocation succeeded */
#define CHECK_ALLOC(ptr) \
    do { \
        if ((ptr) == NULL) { \
            fprintf(stderr, "[ERROR] %s:%d (%s) - Memory allocation failed\n", \
                    __FILE__, __LINE__, __func__); \
            return NUMERICS_ERR_MEMORY; \
        } \
    } while(0)

/** Standardized convergence check - logs warning if not converged */
#define CHECK_CONVERGENCE(error, tol, iter) \
    do { \
        if ((error) > (tol)) { \
            fprintf(stderr, "[WARN] %s:%d (%s) - Not converged: error=%.3e > tol=%.3e at iter=%d\n", \
                    __FILE__, __LINE__, __func__, (double)(error), (double)(tol), (int)(iter)); \
        } \
    } while(0)

/* ===========================================================================
 * Legacy Compatibility
 * =========================================================================== */

/* Map legacy return codes to new error codes */
#define numerics_error_compatible(ret) \
    ((ret) == 0 ? NUMERICS_OK : \
     (ret) == 1 ? NUMERICS_ERR_MAX_ITER : \
     (ret) < 0 ? (NumericsError)(ret) : NUMERICS_ERR_UNKNOWN)

/* Legacy wrapper - converts to new system but maintains compatibility */
#define numerics_error_std(fn_name, msg) \
    do { \
        fprintf(stderr, "[ERROR] %s:%d (%s) - %s: %s\n", \
                __FILE__, __LINE__, __func__, fn_name, msg); \
    } while(0)

#ifdef __cplusplus
}
#endif

#endif /* NUMERICS_ERRORS_H */
