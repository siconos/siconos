#ifndef MLCP_CST_H
#define MLCP_CST_H

enum MLCP_SOLVER {
  SICONOS_MLCP_PGS = 100,
  SICONOS_MLCP_RPGS = 101,
  SICONOS_MLCP_PSOR = 102,
  SICONOS_MLCP_RPSOR = 103,
  SICONOS_MLCP_PATH = 104,
  SICONOS_MLCP_ENUM = 105,
  SICONOS_MLCP_SIMPLEX = 106,
  SICONOS_MLCP_DIRECT_ENUM = 107,
  SICONOS_MLCP_PATH_ENUM = 108,
  SICONOS_MLCP_DIRECT_SIMPLEX = 109,
  SICONOS_MLCP_DIRECT_PATH = 110,
  SICONOS_MLCP_DIRECT_PATH_ENUM = 111,
  SICONOS_MLCP_FB = 112,
  SICONOS_MLCP_DIRECT_FB = 113,
  SICONOS_MLCP_PGS_SBM = 114,
  SICONOS_MLCP_LCP_LEMKE = 115
};

enum SICONOS_IPARAM_MLCP {
  SICONOS_IPARAM_MLCP_PGS_EXPLICIT = 2,
  SICONOS_IPARAM_MLCP_PGS_SUM_ITER = 3,
  SICONOS_IPARAM_MLCP_ENUM_USE_DGELS =
      4,  // activate to use dgels rather than dgesv in mlcp driver (enum only indeed)
  SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS = 5,  // number of possible configurations
  SICONOS_IPARAM_MLCP_UPDATE_REQUIRED = 8,           // true if the problem needs update
};

enum SICONOS_DPARAM_MLCP {
  SICONOS_DPARAM_MLCP_PGS_SUM_ERRORS = 2,
  SICONOS_DPARAM_MLCP_RHO = 3,
  SICONOS_DPARAM_MLCP_OMEGA = 4,
  SICONOS_DPARAM_MLCP_SIGN_TOL_NEG =
      5,  // tolerance for the direct solver, used to check complementarity
  SICONOS_DPARAM_MLCP_SIGN_TOL_POS =
      6,  // tolerance for the direct solver, used to check complementarity
};

extern const char* const SICONOS_NONAME_STR;
extern const char* const SICONOS_MLCP_PGS_STR;
extern const char* const SICONOS_MLCP_RPGS_STR;
extern const char* const SICONOS_MLCP_PSOR_STR;
extern const char* const SICONOS_MLCP_RPSOR_STR;
extern const char* const SICONOS_MLCP_PATH_STR;
extern const char* const SICONOS_MLCP_ENUM_STR;
extern const char* const SICONOS_MLCP_SIMPLEX_STR;
extern const char* const SICONOS_MLCP_DIRECT_ENUM_STR;
extern const char* const SICONOS_MLCP_PATH_ENUM_STR;
extern const char* const SICONOS_MLCP_DIRECT_SIMPLEX_STR;
extern const char* const SICONOS_MLCP_DIRECT_PATH_STR;
extern const char* const SICONOS_MLCP_DIRECT_PATH_ENUM_STR;
extern const char* const SICONOS_MLCP_FB_STR;
extern const char* const SICONOS_MLCP_DIRECT_FB_STR;
extern const char* const SICONOS_MLCP_PGS_SBM_STR;
extern const char* const SICONOS_MLCP_LCP_LEMKE_STR;

#endif
