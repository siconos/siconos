                      Brief Description of
         ODEPACK - A Systematized Collection of ODE Solvers
                    Double Precision Version


                        Alan C. Hindmarsh
          Center for Applied Scientific Computing, L-561
              Lawrence Livermore National Laboratory
                     Livermore, CA 94551, U.S.A.

                           20 June 2001



Work performed under the auspices of the U.S. Department of Energy
by the Lawrence Livermore National Laboratory under contract
No. W-7405-Eng-48, and supported (formerly) by the DOE Office of 
Energy Research, Applied Mathematical Sciences Research Program.

---------------------------------------------------------------------------

   ODEPACK is a collection of Fortran solvers for the initial value
problem for ordinary differential equation systems.  It consists of nine
solvers, namely a basic solver called LSODE and eight variants of it --
LSODES, LSODA, LSODAR, LSODPK, LSODKR, LSODI, LSOIBT, and LSODIS.
The collection is suitable for both stiff and nonstiff systems.  It
includes solvers for systems given in explicit form, dy/dt = f(t,y),
and also solvers for systems given in linearly implicit form, 
A(t,y) dy/dt = g(t,y).  Two of the solvers use general sparse matrix
solvers for the linear systems that arise.  Two others use iterative
(preconditioned Krylov) methods instead of direct methods for these
linear systems.  The most recent addition is LSODIS, which solves
implicit problems with general sparse treatment of all matrices involved.

   The ODEPACK solvers are written in standard Fortran 77, with a few
exceptions, and with minimal machine dependencies.  There are separate
double and single precision versions of ODEPACK.  The actual solver
names are those given above with a prefix of D- or S- for the double
or single precision version, respectively, i.e. DLSODE/SLSODE, etc.
Each solver consists of a main driver subroutine having the same name
as the solver and some number of subordinate routines.  For each
solver, there is also a demonstration program, which solves one or two
simple problems in a somewhat self-checking manner.

   Recently, the ODEPACK solvers were upgraded to improve their
portability in numerous ways.  Among the improvements are (a) renaming
of routines and Common blocks to distinguish double and single
precision versions, (b) use of generic intrinsic function names, (c)
elimination of the Block Data subprogram, (d) use of a portable
routine to set the unit roundoff, and (e) passing of quoted strings to
the error message handler.  In addition, the prologue and internal
comments were reformatted, and use mixed upper/lower case.  Numerous
minor corrections and improvements were also made.  

   The above upgrade operations were applied to LSODE earlier than they
were to the rest of ODEPACK, and the two upgrades were done somewhat
independently.  As a result, some differences will be apparent in the
source files of LSODE and the other solvers -- primarily in the
formatting of the comment line prologue of the main driver routine.
In Subroutines DLSODE/SLSODE and their subordinate routines, the
prologue was written in "SLATEC format", while for the other solvers a
more relaxed style was used.  The differences are entirely cosmetic,
however, and do not affect performance.

   Documentation on the usage of each solver is provided in the
initial block of comment lines in the source file, which (in most
cases) includes a simple example.  A demonstration program (in
seperate double/single precision versions) is also available.

   What follows is a summary of the capabilities of ODEPACK, comments
about usage documentation, and notes about installing the collection.
For additional documentation on ODEPACK, see also the papers [1], [2]
(for LSODE), and [3] (for LSODPK and LSODKR), and in the references
cited there.  (However, the document [2] does not reflect the upgrade
operations described above.)


References:

[1]  A. C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE Solvers,"
     in Scientific Computing, R. S. Stepleman et al. (eds.), North-Holland,
     Amsterdam, 1983 (vol. 1 of IMACS Transactions on Scientific Computation),
     pp. 55-64.

[2]  K. Radhakrishnan and A. C. Hindmarsh, "Description and Use of LSODE,
     the Livermore Solver for Ordinary Differential Equations," LLNL
     report UCRL-ID-113855, December 1993.

[3]  P. N. Brown and A. C. Hindmarsh, "Reduced Storage Matrix Methods
     in Stiff ODE Systems," J. Appl. Math. & Comp., 31 (1989), pp.40-91.

---------------------------------------------------------------------------

                     I. Summary of the ODEPACK Solvers

A. Solvers for explicitly given systems.

For each of the following solvers, it is assumed that the ODEs are
given explicitly, so that the system can be written in the form
dy/dt = f(t,y), where y is the vector of dependent variables, and t is
the independent variable.

1. LSODE (Livermore Solver for Ordinary Differential Equations) is the
   basic solver of the collection.  It solves stiff and nonstiff systems
   of the form dy/dt = f.  In the stiff case, it treats the Jacobian 
   matrix df/dy as either a dense (full) or a banded matrix, and as either 
   user-supplied or internally approximated by difference quotients.  
   It uses Adams methods (predictor-corrector) in the nonstiff case, 
   and Backward Differentiation Formula (BDF) methods (the Gear methods)
   in the stiff case.  The linear systems that arise are solved by direct
   methods (LU factor/solve).  LSODE supersedes the older GEAR and GEARB 
   packages, and reflects a complete redesign of the user interface 
   and internal organization, with some algorithmic improvements.

2. LSODES, written jointly with A. H. Sherman, solves systems dy/dt = f
   and in the stiff case treats the Jacobian matrix in general sparse
   form.  It determines the sparsity structure on its own, or optionally
   accepts this information from the user.  It then uses parts of the
   Yale Sparse Matrix Package (YSMP) to solve the linear systems that
   arise, by a sparse (direct) LU factorization/backsolve method.
   LSODES supersedes, and improves upon, the older GEARS package.

3. LSODA, written jointly with L. R. Petzold, solves systems dy/dt = f
   with a dense or banded Jacobian when the problem is stiff, but it
   automatically selects between nonstiff (Adams) and stiff (BDF) 
   methods.  It uses the nonstiff method initially, and dynamically 
   monitors data in order to decide which method to use.

4. LSODAR, also written jointly with L. R. Petzold, is a variant of 
   LSODA with a rootfinding capability added.  Thus it solves problems 
   dy/dt = f with dense or banded Jacobian and automatic method 
   selection, and at the same time, it finds the roots of any of a 
   set of given functions of the form g(t,y).  This is often useful 
   for finding stop conditions, or for finding points at which a switch
   is to be made in the function f.

5. LSODPK, written jointly with Peter N. Brown, is a variant of LSODE
   in which the direct solvers for the linear systems have been replaced
   by a selection of four preconditioned Krylov (iterative) solvers.  
   The user must supply a pair of routine to evaluate, preprocess, and
   solve the (left and/or right) preconditioner matrices.  LSODPK also
   includes an option for a user-supplied linear system solver to be used
   without Krylov iteration.

6. LSODKR is a variant of LSODPK with the addition of the same
   rootfinding capability as in LSODAR, and also of automatic switching
   between functional and Newton iteration.  The nonlinear iteration 
   method-switching differs from the method-switching in LSODA and LSODAR,
   but provides similar savings by using the cheaper method in the non-stiff
   regions of the problem.  LSODKR also improves on the Krylov methods in
   LSODPK by offering the option to save and reuse the approximate Jacobian
   data underlying the preconditioner.


B. Solvers for linearly implicit systems.

The following solvers treat systems in the linearly implicit form
A(t,y) dy/dt = g(t,y), A = a square matrix, i.e. with the derivative
dy/dt implicit, but linearly so.  These solvers allow A to be
singular, in which case the system is a differential-algebraic
equation (DAE) system.  In that case, the user must be very careful
to supply a well-posed problem with consistent initial conditions.

7. LSODI, written jointly with J. F. Painter, solves linearly implicit
   systems in which the matrices involved (A, dg/dy, and d(A dy/dt)/dy) 
   are all assumed to be either dense or banded.  LSODI supersedes the 
   older GEARIB solver and improves upon it in numerous ways.

8. LSOIBT, written jointly with C. S. Kenney, solves linearly implicit
   systems in which the matrices involved are all assumed to be
   block-tridiagonal.  Linear systems are solved by the LU method.

9. LSODIS, written jointly with S. Balsdon, solves linearly implicit
   systems in which the matrices involved are all assumed to be sparse.
   Like LSODES, LSODIS either determines the sparsity structure or
   accepts it from the user, and uses parts of the Yale Sparse Matrix
   Package to solve the linear systems that arise, by a direct method.

---------------------------------------------------------------------------

                          II. Usage Documentation

   Each of the solvers in the ODEPACK collection is headed by a
user-callable driver subroutine, with the same name as the solver
(DLSODE, etc.).  The call sequence of the driver routine includes the
names of one or more user-supplied subroutines that define the ODE
system, and various other problem and solution parameters.  Complete
user documentation is given in the initial block of comment lines 
(the prologue) in the driver routine.  In each case, this prologue is
organized as follows:

 * Summary of Usage (short, for standard modes of use)
 * Example Problem, with code and output (except for LSODPK and LSODKR)
 * Full Description of User Interface, further divided as follows:
      1. Call sequence description (including optional inputs/outputs)
      2. Optionally callable routines
      3. Descriptions of internal Common blocks
      4. Optionally user-replaceable routines
 * Revision History, showing date written and dates of revisions
 * Other Routines, a list of all subordinate routines for the solver

   First-time users should read only the Summary of Usage and look at the
the Example Problem (or demonstration program), then later refer to the
Full Description if and when more details or nonstandard options are needed.

---------------------------------------------------------------------------

                          III. Installation Notes

1. In addition to this document, the double precision version of ODEPACK 
consists of three source files, plus a demonstration program file.
The solver source files are organized as follows:

   opkdmain.f = Main Source File, consisting of the driver subroutine
                for each of the nine solvers

   opkda1.f   = First Auxiliary Source File, consisting of subordinate
                routines for the solvers

   opkda2.f   = Second Auxiliary Source File, consisting of subordinate
                routines which may already reside on the user's system
                (See Notes 2 and 3 below.)

The demonstration program file is:

   opkddemos  = a merge of the nine demonstration programs, with the
                source for each followed by a sample output

2. The Second Auxiliary Source File includes the routines from the
LINPACK and BLAS collections that are needed by the solvers (and by two
of the demonstration programs), for the solution of dense and banded
linear systems and associated basic linear algebra operations.
These routine are:
   From LINPACK:  DGEFA, DGESL, DGBFA, DGBSL
   From the BLAS: DAXPY, DCOPY, DDOT, DSCAL, DNRM2, IDAMAX
If your computer system already has these routines, and especially if it
has machine-optimized versions, the copies provided here can be discarded.

3. The Second Auxiliary Source File includes a set of five routines --
XERRWD, XSETUN, XSETF, IXSAV, IUMACH -- which handle error messages
from the solvers.  This set is in fact a reduced version (sufficient
for the needs of ODEPACK) of a much larger error handling package from
the SLATEC Library.  If your computer system already has the full
SLATEC error handler, the version provided here can be discarded.  If
the reduced version is used, its machine-dependent features should be
checked first; see comments in Subroutine XERRWD.

4. ODEPACK contains a few instances where ANSI Fortran 77 is violated:
   (a) In various places in the LSODES and LSODIS solvers, a call to a
       subroutine has a subscripted real array as an argument where the
       subroutine called expects an integer array.  Calls of this form
       occur in Subroutine DLSODES (to DSTODE), in DIPREP (to DPREP),
       in Subroutine DLSODIS (to DSTODI), and in DIPREPI (to DPREPI).
       Another such call occurs in the DLSODES demonstration program, 
       from the main program to Subroutine SSOUT.  This is done in order
       to use work space in an efficient manner, as the same space is
       sometimes used for real work space and sometimes for integer work
       space.  If your compiler does not accept this feature, one possible
       way to get the desired result is to compile the called routines
       and calling routines in separate jobs, and then combine the binary
       modules in an appropriate manner.  If this procedure is still not
       acceptable under your system, it will be necessary to radically
       alter the structure of the array RWORK within the LSODES or LSODIS
       solver package.  (See also Note 5 below.)
   (b) Each ODEPACK solver treats the arguments NEQ, Y, RTOL, and ATOL 
       as arrays, even though the length may be only 1.  Moreover, 
       except for Y, the usage instructions say that these arguments 
       may be either arrays or scalars.  If your system does not allow 
       such a mismatch, then the documentation of these arguments 
       should be changed accordingly.

5. For maximum storage economy, the LSODES and LSODIS solvers make use
of the real to integer wordlength ratio.  This is assumed to be an
integer L such that if a real array R and an integer array M occupy
the same space in memory, R(1) having the same bit address as M(1),
then R(I) has the same address as M((I-1)*L+1).  This ratio L is
usually 2 for double precision, and this is the value used in the
double precision version supplied.  If this value is incorrect, it
needs to be changed in two places:
  (a) The integer LENRAT is DATA-loaded in Subroutines DLSODES and
      DLSODIS to this ratio, shortly below the prologue.
  (b) The integer LRATIO is DATA-loaded in Subroutine CDRV to this
      ratio, shortly below the prologue of that routine.
(See comments in both places.)  If the ratio is not an integer, use 
the greatest integer not exceeding the ratio.

6. For installation of ODEPACK on a Cray computer, the source files
supplied include compiler directives for the CFT compiler.  These have
the form CDIR$ IVDEP and occur prior to certain loops that involve
subscript shifts (and would otherwise not be vectorized).  These
directives are (or should be) treated as comments by any other compiler.

7. On first obtaining ODEPACK, the demonstration programs should be
compiled and executed prior to any other use of the solvers.  In most
cases, these excercise all of the major method options in each solver,
and are self-checking.  (In the case of LSODPK and LSODKR, the
demonstration programs are not self-checking, and for LSODKR only one
major method option is used.)  In any case, the output can be compared
with the sample output supplied, which was generated from the double
precision version of ODEPACK on a 32-bit computer.  When comparing your
output with that supplied, differences of 10-20% in the final values of
the various statistical counters can be attributed to differences in
the roundoff properties of different computer systems.

8. If some subset of the whole ODEPACK collection is desired, without
unneeded routines, the appropriate routines must be extracted
accordingly.  The following lists give the routines needed for the
double precision version of each solver.

    The DLSODE solver consists of the routines
  DLSODE, DINTDY, DSTODE, DCFODE, DPREPJ, DSOLSY, DEWSET, DVNORM, DSRCOM,
  DGEFA, DGESL, DGBFA, DGBSL, DAXPY, DSCAL, DDOT, IDAMAX,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH

    The DLSODES solver consists of the routines
  DLSODES, DIPREP, DPREP, JGROUP, ADJLR, CNTNZU, DINTDY, DSTODE, DCFODE,
  DPRJS, DSOLSS, DEWSET, DVNORM, DSRCMS,
  ODRV, MD, MDI, MDM, MDP, MDU, SRO,
  CDRV, NROC, NSFC, NNFC, NNSC, NNTC,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH

    The DLSODA solver consists of the routines
  DLSODA, DINTDY, DSTODA, DCFODE, DPRJA, DSOLSY, DEWSET,
  DMNORM, DFNORM, DBNORM, DSRCMA,
  DGEFA, DGESL, DGBFA, DGBSL, DAXPY, DSCAL, DDOT, IDAMAX,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODAR solver consists of the routines
  DLSODAR, DRCHEK, DROOTS, DINTDY, DSTODA, DCFODE, DPRJA, DSOLSY, DEWSET,
  DMNORM, DFNORM, DBNORM, DSRCAR,
  DGEFA, DGESL, DGBFA, DGBSL, DAXPY, DSCAL, DDOT, DCOPY, IDAMAX,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODPK solver consists of the routines
  DLSODPK, DINTDY, DEWSET, DVNORM, DSTODPK, DCFODE, DPKSET, DSOLPK,
  DSPIOM, DATV, DORTHOG, DHEFA, DHESL, DSPIGMR, DHEQR, DHELS, 
  DPCG, DPCGS, DATP, DUSOL, DSRCPK,
  DAXPY, DSCAL, DCOPY, DDOT, DNRM2, IDAMAX,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODKR solver consists of the routines
  DLSODKR, DRCHEK, DROOTS, DLHIN, DINTDY, DEWSET, DVNORM, DSTOKA,
  DCFODE, DSETPK, DSOLPK, DSPIOM, DATV, DORTHOG, DHEFA, DHESL, DSPIGMR,
  DHEQR, DHELS, DPCG, DPCGS, DATP, DUSOL, DSRCKR,
  DAXPY, DSCAL, DCOPY, DDOT, DNRM2, IDAMAX,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSODI solver consists of the routines
  DLSODI, DAINVG, DINTDY, DSTODI, DCFODE, DPREPJI, DSOLSY, DEWSET,
  DVNORM, DSRCOM,
  DGEFA, DGESL, DGBFA, DGBSL, DAXPY, DSCAL, DDOT, IDAMAX,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH

     The DLSOIBT solver consists of the routines
  DLSOIBT, DAIGBT, DINTDY, DSTODI, DCFODE, DPJIBT, DSLSBT, DEWSET,
  DVNORM, DSRCOM, DDECBT, DSOLBT,
  DGEFA, DGESL, DAXPY, DSCAL, DDOT, IDAMAX,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH
 
     The DLSODIS solver consists of the routines
  DLSODIS, DAINVGS, DIPREPI, DPREPI, JGROUP, ADJLR, CNTNZU, DINTDY,
  DSTODI, DCFODE, DPRJIS, DSOLSS, DEWSET, DVNORM, DSRCMS,
  ODRV, MD, MDI, MDM, MDP, MDU, SRO,
  CDRV, NROC, NSFC, NNFC, NNSC, NNTC,
  DUMACH, XERRWD, XSETUN, XSETF, IXSAV, IUMACH
