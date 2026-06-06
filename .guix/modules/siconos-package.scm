;;; SPDX-License-Identifier: GPL-3.0-or-later

(define-module (siconos-package)
  #:use-module (guix)
  #:use-module (guix build-system cmake)
  #:use-module (guix build-system copy)
  #:use-module (guix build-system pyproject)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix packages)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages boost)
  #:use-module (gnu packages check)
  #:use-module (gnu packages cmake)
  #:use-module (gnu packages commencement)
  #:use-module (gnu packages cpp)
  #:use-module ((gnu packages game-development) #:select (bullet))
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages llvm)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages mpi)
  #:use-module (gnu packages multiprecision)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages pretty-print)
  #:use-module (gnu packages python)
  #:use-module (gnu packages python-build)
  #:use-module (gnu packages python-science)
  #:use-module (gnu packages python-xyz)
  #:use-module (gnu packages swig)
  #:use-module (gnu packages version-control)
  #:use-module (gnu packages xml)
  #:use-module (guix-science packages physics)
  #:use-module (guix-science-nonfree packages cuda))


;(define vcs-file?
;  ;; Renvoie vrai lorsque le fichier donné est sous contrôle de version.
;  (or (git-predicate (dirname (dirname (current-source-directory))))
;      (const #t)))                                ; pas dans un dépôt Git
(define-public fclib-devel
  (package
    (name "fclib-devel")
    (version "3.1.0")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/FrictionalContactLibrary/fclib/")
             (commit "f894f34f02093910ceb49943d69f4901142bb880")))
       (sha256
        (base32 "0dg2r54jl8sr8pkc0r0q5sx0bwm0z093dkfcfci67x51z5w8nz98"))))
    (build-system cmake-build-system)
    (arguments
     (list
      #:configure-flags
      #~(list "-DFCLIB_HEADER_ONLY=OFF")))
    (propagated-inputs (list hdf5))
    (home-page "https://frictionalcontactlibrary.github.io/")
    (synopsis "A collection of discrete 3D Frictional Contact (FC) problems")
    (description
     "FCLIB is an open source collection of Frictional
Contact (FC) problems stored in a specific HDF5 format with a light
implementation in C Language of Input/Output functions to read and write those
problems.")
    (license license:asl2.0)))


(define-public siconos-devel
  (package
    (name "siconos-devel")
    (version "x")
    (source
     (local-file "../.." "siconos-checkout"
                 #:recursive? #t
     ;            #:select? vcs-file?
     ))
    
    (build-system cmake-build-system)
    (arguments
     (list
        #:tests? #f
	#:modules `(
		    (guix build cmake-build-system)
		    ((guix build pyproject-build-system) #:prefix py:)
		    (guix build utils)
		   )
        #:imported-modules (append %cmake-build-system-modules
				   %pyproject-build-system-modules)
        #:configure-flags
        #~(list "-DCMAKE_VERBOSE_MAKEFILE=1"
                "-DWITH_GUIX=1"
                (string-append "-Dpetsc_ROOT="
                               #$(this-package-input "petsc-openmpi"))
                (string-append "-DHDF5_ROOT="
                               #$(this-package-input "hdf5-parallel-openmpi") "/lib")
                (string-append "-DSICONOS_CUSTOM_INSTALL=" #$output)
                "-DCOMPONENTS=externals;numerics;kernel;control;mechanics;io;vkernel"
                "-DCMAKE_BUILD_TYPE=Release"
                "-DWITH_PYB11_WRAPPER=1"
                "-DFCLIB_ROOT=1"
                "-DWITH_BULLET=1"
                "-DWITH_OCC=1"
                "-DWITH_FCLIB=1"
                "-DISOLATED_INSTALL=1"
                "-DWITH_TESTING=0"
                "-DWITH_vkernel_TESTING=0"
                "-DWITH_BULLET=1"
                "-DWITH_HDF5=1"
                "-DWITH_PETSC=1"
                "-DWITH_MPI=1"
                "-DWITH_OPENMP=1")
        #:phases
        #~(modify-phases %standard-phases
            (add-after 'check 'set-SOURCE-DATE-EPOCH
              (lambda _
                (setenv "SOURCE_DATE_EPOCH" "315532800")))
            (add-after 'unpack 'some-quick-patches
              (lambda _
                ;; compiler specification not needed.
                (substitute* "siconos-config.cmake.in"
                  (("set\\(CMAKE_CXX_COMPILER @CMAKE_CXX_COMPILER@\\)")
                   "")
                  (("set\\(CMAKE_C_COMPILER @CMAKE_C_COMPILER@\\)")
                   "")
                  (("set\\(CMAKE_Fortran_COMPILER @CMAKE_Fortran_COMPILER@\\)")
                   ""))

                ;; failure with recent compilers
                (substitute* "externals/numeric_bindings/boost/numeric/bindings/blas/detail/cblas.h"
                  (("#ifdef HAS_OpenBLAS")
                   "#ifdef REMOVED_HAS_OpenBLAS"))

                ;; guix needs CMAKE_INSTALL_PREFIX and it is ok.
                (substitute* "cmake/SiconosSetup.cmake"
                  (("FATAL_ERROR")
                   "WARNING"))

                ;; other non-fatal fatal error
                (substitute* "cmake/SiconosInstallSetup.cmake"
                  (("FATAL_ERROR")
                   "WARNING"))

                ;; Some of the Siconos/numerics tests are unstable.
                (substitute* "numerics/numerics_tests.cmake"
                  (("WITH_TESTING")
                   "WITH_UNSTABLE_TESTING"))

                ;; Cable is a new feature and is not yet stable
                (substitute* "mechanics/mechanics_tests.cmake"
                  (("begin_tests\\(src/fem/cable/test")
                   "# begin_tests(src/fem/cable/test")
                  (("new_test\\(SOURCES CableDSTest.cpp")
                   "# new_test(SOURCES CableDSTest.cpp"))))

                ;; wrap siconos file
            	    (add-after 'install 'wrap-python
              (assoc-ref py:%standard-phases 'wrap))
	      	    )))

    (native-inputs (list clang-toolchain-21
                         pkg-config
                         cppunit
                         gfortran
                         git-minimal/pinned
                         pybind11
                         python-lxml
                         python-pytest
                         swig))
      (inputs (list python))
      (propagated-inputs (list boost
                               bullet
                               eigen
                               fclib-devel
                               gmp
                               hdf5-parallel-openmpi
                               lapack
                               nlohmann-json
                               openblas
                               opencascade-occt
                               petsc-openmpi
                               openmpi-4
                               python-h5py
                               python-numpy
                               python-packaging
                               python-petsc4py
                               python-scipy
                               python-occ-core
                               python-wheel
                               suitesparse))
      (home-page
       "https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/index.html")
      (synopsis "Library for nonsmooth numerical simulation")
      (description
        "Siconos is an open-source scientific software primarily targeted at
modeling and simulating nonsmooth dynamical systems in C++ and in Python:
Mechanical systems (rigid or solid) with unilateral contact and Coulomb
friction and impact (nonsmooth mechanics, contact dynamics, multibody systems
dynamics or granular materials).  Switched Electrical Circuit such as
electrical circuits with ideal and piecewise linear components: power
converter, rectifier, Phase-Locked Loop (PLL) or Analog-to-Digital converter.
Sliding mode control systems.  Biology (Gene regulatory network).

Other applications are found in Systems and Control (hybrid systems,
differential inclusions, optimal control with state constraints),
Optimization (Complementarity systems and Variational inequalities), Fluid
Mechanics, and Computer Graphics.")
      (license license:asl2.0)))

(define-public siconos-devel-cuda
  (package
    (inherit siconos-devel)
    (name "siconos-devel-cuda")
    (native-inputs (modify-inputs (package-native-inputs siconos-devel)
                     ;; prepend gcc-toolchain-14 to ensure it takes
                     ;; precedence for cuda compilation, as cuda
                     ;; requires matching host compiler and standard
                     ;; library versions, with append it fails here.
                     (prepend gcc-toolchain-14)))
    (inputs (modify-inputs (package-inputs siconos-devel)
              (prepend cuda)))
    (arguments
     (substitute-keyword-arguments (package-arguments siconos-devel)
       ((#:configure-flags flags '())
        #~(append (list
                   "-DCMAKE_CXX_COMPILER=clang++"
                   "-DWITH_GPU=1"
                   "-DWITH_CUDA=1"
                   "-DCMAKE_CUDA_ARCHITECTURES=80;86;90"
                   "-DCMAKE_CUDA_HOST_COMPILER=gcc")
                  #$flags))))))

siconos-devel
