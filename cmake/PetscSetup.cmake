# Find petsc using pkg-config
if(DEFINED petsc_ROOT)
  set(ENV{PKG_CONFIG_PATH} ${petsc_ROOT}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH})
elseif(DEFINED ENV{petsc_ROOT})
  set(ENV{PKG_CONFIG_PATH} $ENV{petsc_ROOT}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH})
endif()

include(FindPkgConfig)

set(PKG_CONFIG_ARGN "--define-variable=PETSC_DIR=/tmp")

pkg_search_module(PETSC
  REQUIRED
  IMPORTED_TARGET
  petsc PETSc)

#find_package(MPI REQUIRED)
message("- includes for petsc ${PETSC_INCLUDEDIR}")
message("- lib path for petsc ${PETSC_LIBDIR}")
message("- libraries for petsc ${PETSC_LINK_LIBRARIES}")
message("- flags for petsc ${PETSC_CFLAGS}")
message("- ldflags for petsc ${PETSC_LDFLAGS}")

if(WITH_HDF5)
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 COMPONENTS C REQUIRED)
endif()

