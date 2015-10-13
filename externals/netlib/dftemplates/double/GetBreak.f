*     ================================================================
      DOUBLE PRECISION FUNCTION GETBREAK()
*
*     Get breakdown parameter tolerance; for the test routine,
*     set to machine precision.
*
      DOUBLE PRECISION EPS, DLAMCH
*
      EPS = DLAMCH('EPS')
      GETBREAK = EPS**2
*
      RETURN
*
      END
