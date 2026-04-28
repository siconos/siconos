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
  #:use-module (gnu packages cpp)
  #:use-module ((gnu packages game-development) #:select (bullet))
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages llvm)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages multiprecision)
  #:use-module (gnu packages pretty-print)
  #:use-module (gnu packages python)
  #:use-module (gnu packages python-build)
  #:use-module (gnu packages python-science)
  #:use-module (gnu packages python-xyz)
  #:use-module (gnu packages swig)
  #:use-module (gnu packages version-control)
  #:use-module (gnu packages xml)
  #:use-module (guix-science packages physics)
  )


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
        #:imported-modules `((guix build python-build-system)
                             ,@%cmake-build-system-modules)
        #:configure-flags
        #~(list "-DCMAKE_VERBOSE_MAKEFILE=1"
                "-DWITH_PYB11_WRAPPER=1"
                "-DFCLIB_ROOT=1"
                "-DWITH_BULLET=1"
                "-DWITH_OCC=1"
                "-DWITH_FCLIB=1"
                "-DISOLATED_INSTALL=1"
                "-DWITH_TESTING=0"
                "-DWITH_vkernel_TESTING=0"
                (string-append "-DSICONOS_CUSTOM_INSTALL=" #$output)
                "-DCOMPONENTS=externals;numerics;kernel;control;mechanics;io;vkernel")
        #:phases
        #~(modify-phases %standard-phases
            (add-after 'check 'set-SOURCE-DATE-EPOCH
              (lambda _
                (setenv "SOURCE_DATE_EPOCH" "315532800")))
            (add-after 'unpack 'some-quick-patches
              (lambda _
                (substitute* "cmake/fclib_setup.cmake"
                  ;; compiler specfication fails under guix, but it is not needed.
                  (("target_compile_options\\(fclib INTERFACE")
                   "message(STATUS \"pffff.\" ")
                  (("\\$<\\$<CXX_COMPILER_ID")
                   "#")
                  ;; fclib config not found.
                  (("find_package\\(FCLIB 3.0.0 CONFIG REQUIRED\\)")
                   "find_package(FCLIB 3.0.0 CONFIG REQUIRED)
set(ConfigPackageLocation lib/cmake/siconos-${SICONOS_VERSION})"))
                
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

            ;; installation of python files with setup.py fails with
            ;; recents python/numpy
            (add-after 'install 'install-python-files
              (lambda* (#:key inputs outputs #:allow-other-keys)
                (let* ((python-version (@ (guix build python-build-system)
                                          python-version))
                       (out (assoc-ref outputs "out"))
                       (version (python-version (assoc-ref inputs "python")))
                       (pydir (string-append out "/lib/python" version
                                             "/site-packages/")))
                  (chdir "./python")
                  (for-each (lambda (file)
                              (let ((dirinst (string-append pydir "/"
                                                            (dirname file))))
                                (mkdir-p dirinst)
                                (install-file file dirinst)))
                            (find-files "siconos" "\\.py$"))))))))

      (native-inputs (list clang-toolchain-21
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
                               fmt
                               gmp
                               lapack
                               nlohmann-json
                               openblas
                               opencascade-occt
                               python-h5py
                               python-numpy
                               python-packaging
                               python-scipy
                               python-occ-core
                               python-wheel
                               range-v3
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

siconos-devel
