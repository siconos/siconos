

#define WITH_CMAKE

#ifndef SVN_REVISION
#cmakedefine SVN_REVISION ${SVN_REVISION}
#endif

#cmakedefine HAVE_BLAS
#cmakedefine HAVE_LAPACK
#cmakedefine HAVE_ATLAS
#cmakedefine HAVE_CLAPACK_H
#cmakedefine HAVE_CBLAS_H
#cmakedefine HAVE_XERBLA
#cmakedefine HAVE_PATHFERRIS
#cmakedefine HAVE_MLCPSIMPLEX
