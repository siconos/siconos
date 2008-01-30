PROJECT(dftemplates)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(dftemplates_library_SRCS    
                        Utils.f 
			MatVec.f PSolve.f 
			CG.f CGREVCOM.f 
			BiCG.f BiCGREVCOM.f 
			BiCGSTAB.f BiCGSTABREVCOM.f 
			CGS.f CGSREVCOM.f 
			Cheby.f ChebyREVCOM.f 
			GMRES.f GMRESREVCOM.f 
			Jacobi.f JacobiREVCOM.f 
			QMR.f QMRREVCOM.f 
			SOR.f SORREVCOM.f 
			STOPTEST2.f
			GetBreak.f Matnorm.f)

SET(dftemplates_HDRS dftemplates.h)

#F95 Bug
CONFIGURE_FILE(dftemplates.h saved.mod.proxy COPYONLY)
CONFIGURE_FILE(dftemplates.h dgemv..mod.proxy COPYONLY)

# actual target:
ADD_LIBRARY(dftemplates ${dftemplates_library_HDRS} ${dftemplates_library_SRCS})
SET_TARGET_PROPERTIES(dftemplates PROPERTIES VERSION 0.0.0)

# add install target:
INSTALL(TARGETS dftemplates DESTINATION lib)
