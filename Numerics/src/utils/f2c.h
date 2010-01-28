/* Siconos-Numerics, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
/* This is also plain -*-c-*-
// Copyright (C) 2004
// Christian Stimming <stimming@tuhh.de>

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2, or (at
// your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License along
// with this library; see the file COPYING.  If not, write to the Free
// Software Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307,
// USA.
*/

/*
// NOTE: this is a modified file of f2c.h to work with most C++ compilers.
//       f2c.h for example, defines abs() as a macro, casuing parsing
//       problems when it encounters abs() in stdlib.h, and so on.
//       It is needed in Lapack++ because some of the C/Lapack functions, like
//       ilaenv_(), need to know what ftn_len is (system dependent.)
//
//       pow_ri() has been removed, since it collides with
//       sunmath.h on SUN4s.
//
//        "typedef float real" has been removed, due to conflict with
//        AT&T CC with real() for complex.
//
//        "typedef ... complex" has been renamed to f2c_complex to
//        avoid conflicts with C++ complex class.
*/

/** @file utils/f2c.h
 * @brief Standard Fortran to C header file
 *
 * This file includes many declarations of types and headers that are
 * needed for FORTRAN compatibility of the C code. */

#ifndef _F2C_INCLUDE_H
#define _F2C_INCLUDE_H

typedef long int integer;
typedef char *address;
typedef short int shortint;
typedef float f2c_real;
typedef double doublereal;
typedef struct
{
  f2c_real r, i;
} f2c_complex;

/*
The complex type inside LAPACK++ has to be compatible to the FORTRAN
type, since LAPACK++ has to pass these values back and forth to the
FORTRAN functions. Don't use this data type anywhere outside
anymore. The included C++ class la::complex provides
automatic conversion from and to this FORTRAN type. */
typedef struct
{
  doublereal r, i;
} doublecomplex;
typedef long int logical;
typedef short int shortlogical;

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

extern int ungetc();    /* This is not declared in Sun's <stdio.h> */

/* I/O stuff */
#ifndef DOXYGEN_IGNORE

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long flag;
typedef long ftnlen;
typedef long ftnint;
#endif /* f2c_i2 */

/*external read, write*/
typedef struct
{
  flag cierr;
  ftnint ciunit;
  flag ciend;
  char *cifmt;
  ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{
  flag icierr;
  char *iciunit;
  flag iciend;
  char *icifmt;
  ftnint icirlen;
  ftnint icirnum;
} icilist;

/*open*/
typedef struct
{
  flag oerr;
  ftnint ounit;
  char *ofnm;
  ftnlen ofnmlen;
  char *osta;
  char *oacc;
  char *ofm;
  ftnint orl;
  char *oblnk;
} olist;

/*close*/
typedef struct
{
  flag cerr;
  ftnint cunit;
  char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{
  flag aerr;
  ftnint aunit;
} alist;

/* inquire */
typedef struct
{
  flag inerr;
  ftnint inunit;
  char *infile;
  ftnlen infilen;
  ftnint  *inex;  /*parameters in standard's order*/
  ftnint  *inopen;
  ftnint  *innum;
  ftnint  *innamed;
  char    *inname;
  ftnlen  innamlen;
  char    *inacc;
  ftnlen  inacclen;
  char    *inseq;
  ftnlen  inseqlen;
  char    *indir;
  ftnlen  indirlen;
  char    *infmt;
  ftnlen  infmtlen;
  char    *inform;
  ftnint  informlen;
  char    *inunf;
  ftnlen  inunflen;
  ftnint  *inrecl;
  ftnint  *innrec;
  char    *inblank;
  ftnlen  inblanklen;
} inlist;

#define VOID void

union Multitype     /* for multiple entry points */
{
  shortint h;
  integer i;
  f2c_real r;
  doublereal d;
  f2c_complex c;
  doublecomplex z;
};

typedef union Multitype Multitype;

typedef long Long;  /* No longer used; formerly in Namelist */

struct Vardesc      /* for Namelist */
{
  char *name;
  char *addr;
  ftnlen *dims;
  int  type;
};
typedef struct Vardesc Vardesc;

struct Namelist
{
  char *name;
  Vardesc **vars;
  int nvars;
};
typedef struct Namelist Namelist;

/* procedure parameter types for -A and -C++ */
#endif /* DOXYGEN_IGNORE */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint(*J_fp)(...);
typedef integer(*I_fp)(...);
typedef f2c_real(*R_fp)(...);
typedef doublereal(*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical(*L_fp)(...);
typedef shortlogical(*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else /* __cplusplus */
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint(*J_fp)();
typedef integer(*I_fp)();
typedef f2c_real(*R_fp)();
typedef doublereal(*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical(*L_fp)();
typedef shortlogical(*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif /* __cplusplus */
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;   /* (single precision ) complex function */
typedef VOID H_f;   /* character function */
typedef VOID Z_f;   /* double complex function */
typedef doublereal E_f; /* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif /* Skip_f2c_Undefs */


/* If you are using a C++ compiler, append the following to f2c.h
   for compiling libF77 and libI77. */

#ifdef __cplusplus
extern "C"
{
  extern int abort_(void);
  extern double c_abs(f2c_complex *);
  extern void c_cos(f2c_complex *, f2c_complex *);
  extern void c_div(f2c_complex *, f2c_complex *, f2c_complex *);
  extern void c_exp(f2c_complex *, f2c_complex *);
  extern void c_log(f2c_complex *, f2c_complex *);
  extern void c_sin(f2c_complex *, f2c_complex *);
  extern void c_sqrt(f2c_complex *, f2c_complex *);
  extern double d_abs(double *);
  extern double d_acos(double *);
  extern double d_asin(double *);
  extern double d_atan(double *);
  extern double d_atn2(double *, double *);
  extern void d_cnjg(doublecomplex *, doublecomplex *);
  extern double d_cos(double *);
  extern double d_cosh(double *);
  extern double d_dim(double *, double *);
  extern double d_exp(double *);
  extern double d_imag(doublecomplex *);
  extern double d_int(double *);
  extern double d_lg10(double *);
  extern double d_log(double *);
  extern double d_mod(double *, double *);
  extern double d_nint(double *);
  extern double d_prod(float *, float *);
  extern double d_sign(double *, double *);
  extern double d_sin(double *);
  extern double d_sinh(double *);
  extern double d_sqrt(double *);
  extern double d_tan(double *);
  extern double d_tanh(double *);
  extern double derf_(double *);
  extern double derfc_(double *);
  extern integer do_fio(ftnint *, char *, ftnlen);
  extern integer do_lio(ftnint *, ftnint *, char *, ftnlen);
  extern integer do_uio(ftnint *, char *, ftnlen);
  extern integer e_rdfe(void);
  extern integer e_rdue(void);
  extern integer e_rsfe(void);
  extern integer e_rsfi(void);
  extern integer e_rsle(void);
  extern integer e_rsli(void);
  extern integer e_rsue(void);
  extern integer e_wdfe(void);
  extern integer e_wdue(void);
  extern integer e_wsfe(void);
  extern integer e_wsfi(void);
  extern integer e_wsle(void);
  extern integer e_wsli(void);
  extern integer e_wsue(void);
  extern int ef1asc_(ftnint *, ftnlen *, ftnint *, ftnlen *);
  extern integer ef1cmc_(ftnint *, ftnlen *, ftnint *, ftnlen *);
  /*
    #if defined(OS_WIN32) || LAPACK_OS_WIN32 || defined(OS_freebsd)
    extern double erf(double);
    extern double erfc(double);
    #else
    extern double erf(double) throw ();
    extern double erfc(double) throw ();
    #endif
  */
  /* 2004-09-08 Christian Stimming <stimming@tuhh.de>: Remove
     declaration of erf() and erfc() since they cause trouble with
     FreeBSD and they seem to be unneeded anyway. */
  extern double erf_(float *);
  extern double erfc_(float *);
  extern integer f_back(alist *);
  extern integer f_clos(cllist *);
  extern integer f_end(alist *);
  extern void f_exit(void);
  extern integer f_inqu(inlist *);
  extern integer f_open(olist *);
  extern integer f_rew(alist *);
  extern int flush_(void);
  extern void getarg_(integer *, char *, ftnlen);
  extern void getenv_(char *, char *, ftnlen, ftnlen);
  extern short h_abs(short *);
  extern short h_dim(short *, short *);
  extern short h_dnnt(double *);
  extern short h_indx(char *, char *, ftnlen, ftnlen);
  extern short h_len(char *, ftnlen);
  extern short h_mod(short *, short *);
  extern short h_nint(float *);
  extern short h_sign(short *, short *);
  extern short hl_ge(char *, char *, ftnlen, ftnlen);
  extern short hl_gt(char *, char *, ftnlen, ftnlen);
  extern short hl_le(char *, char *, ftnlen, ftnlen);
  extern short hl_lt(char *, char *, ftnlen, ftnlen);
  extern integer i_abs(integer *);
  extern integer i_dim(integer *, integer *);
  extern integer i_dnnt(double *);
  extern integer i_indx(char *, char *, ftnlen, ftnlen);
  extern integer i_len(char *, ftnlen);
  extern integer i_mod(integer *, integer *);
  extern integer i_nint(float *);
  extern integer i_sign(integer *, integer *);
  extern integer iargc_(void);
  extern ftnlen l_ge(char *, char *, ftnlen, ftnlen);
  extern ftnlen l_gt(char *, char *, ftnlen, ftnlen);
  extern ftnlen l_le(char *, char *, ftnlen, ftnlen);
  extern ftnlen l_lt(char *, char *, ftnlen, ftnlen);
  extern void pow_ci(f2c_complex *, f2c_complex *, integer *);
  extern double pow_dd(double *, double *);
  extern double pow_di(double *, integer *);
  extern short pow_hh(short *, shortint *);
  extern integer pow_ii(integer *, integer *);
  /* extern double pow_ri(float *, integer *); */
  extern void pow_zi(doublecomplex *, doublecomplex *, integer *);
  extern void pow_zz(doublecomplex *, doublecomplex *, doublecomplex *);
  extern double r_abs(float *);
  extern double r_acos(float *);
  extern double r_asin(float *);
  extern double r_atan(float *);
  extern double r_atn2(float *, float *);
  extern void r_cnjg(f2c_complex *, f2c_complex *);
  extern double r_cos(float *);
  extern double r_cosh(float *);
  extern double r_dim(float *, float *);
  extern double r_exp(float *);
  extern double r_imag(f2c_complex *);
  extern double r_int(float *);
  extern double r_lg10(float *);
  extern double r_log(float *);
  extern double r_mod(float *, float *);
  extern double r_nint(float *);
  extern double r_sign(float *, float *);
  extern double r_sin(float *);
  extern double r_sinh(float *);
  extern double r_sqrt(float *);
  extern double r_tan(float *);
  extern double r_tanh(float *);
  extern void s_cat(char *, char **, integer *, integer *, ftnlen);
  extern integer s_cmp(char *, char *, ftnlen, ftnlen);
  extern void s_copy(char *, char *, ftnlen, ftnlen);
  extern int s_paus(char *, ftnlen);
  extern integer s_rdfe(cilist *);
  extern integer s_rdue(cilist *);
  extern integer s_rnge(char *, integer, char *, integer);
  extern integer s_rsfe(cilist *);
  extern integer s_rsfi(icilist *);
  extern integer s_rsle(cilist *);
  extern integer s_rsli(icilist *);
  extern integer s_rsne(cilist *);
  extern integer s_rsni(icilist *);
  extern integer s_rsue(cilist *);
  extern int s_stop(char *, ftnlen);
  extern integer s_wdfe(cilist *);
  extern integer s_wdue(cilist *);
  extern integer s_wsfe(cilist *);
  extern integer s_wsfi(icilist *);
  extern integer s_wsle(cilist *);
  extern integer s_wsli(icilist *);
  extern integer s_wsne(cilist *);
  extern integer s_wsni(icilist *);
  extern integer s_wsue(cilist *);
  extern void sig_die(char *, int);
  extern integer signal_(integer *, void (*)(int));
  extern int system_(char *, ftnlen);
  extern double z_abs(doublecomplex *);
  extern void z_cos(doublecomplex *, doublecomplex *);
  extern void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
  extern void z_exp(doublecomplex *, doublecomplex *);
  extern void z_log(doublecomplex *, doublecomplex *);
  extern void z_sin(doublecomplex *, doublecomplex *);
  extern void z_sqrt(doublecomplex *, doublecomplex *);
}
#endif /* __cplusplus  */

#endif /* _F2C_INCLUDE_H */
