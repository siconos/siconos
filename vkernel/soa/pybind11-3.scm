;; https://codeberg.org/guix/guix/pulls/4696

(use-modules ((guix licenses) #:prefix license:)
             (guix packages)
             (guix git-download)
             (guix build-system cmake)
             (guix build-system pyproject)
             (gnu packages cmake)
             (gnu packages python-build)
             (gnu packages python))

(define pybind11-3
  (package
    (name "pybind11")
    (version "3.0.1")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/pybind/pybind11")
             (commit (string-append "v" version))))
       (file-name (git-file-name name version))
       (sha256
        (base32 "1gpax61ndhbr1r179bbgavh50j3aylkddzvbklhyj51mq4d0sb36"))))
    (build-system pyproject-build-system)
    (arguments
     (list
      #:tests? #f  ; Tests require building C++ test modules
      #:phases
      #~(modify-phases %standard-phases
          (add-after 'unpack 'prepare-build
            ;; Upstream removed setup.py in 3.x and now requires scikit-build-core,
            ;; which would create a circular dependency.  Instead, we run CMake
            ;; to populate pybind11/include and pybind11/share, then use a simple
            ;; setup.py.  Based on pybind11 2.13.6's setup.py approach.
            (lambda _
              (delete-file "pyproject.toml")
              ;; Run CMake to install headers and cmake files into pybind11/
              (invoke "cmake" "-S" "." "-B" "build"
                      "-DCMAKE_INSTALL_PREFIX=pybind11"
                      "-DBUILD_TESTING=OFF"
                      "-DPYBIND11_NOPYTHON=ON"
                      "-Dprefix_for_pc_file=${pcfiledir}/../../")
              (invoke "cmake" "--install" "build")
              ;; Write _version.py
              (call-with-output-file "pybind11/_version.py"
                (lambda (port)
                  (format port "__version__ = \"~a\"\n" #$version)
                  (format port "version_info = (~a)\n"
                          (string-join (string-split #$version #\.) ", "))))
              ;; Write setup.py - must explicitly list packages containing
              ;; non-Python files since find_packages() only finds dirs with .py files
              (call-with-output-file "setup.py"
                (lambda (port)
                  (format port "
from setuptools import setup
setup(
    name='pybind11',
    version='~a',
    packages=[
        'pybind11',
        'pybind11.include.pybind11',
        'pybind11.include.pybind11.detail',
        'pybind11.include.pybind11.eigen',
        'pybind11.include.pybind11.stl',
        'pybind11.include.pybind11.conduit',
        'pybind11.share.cmake.pybind11',
        'pybind11.share.pkgconfig',
    ],
    package_data={
        'pybind11': ['py.typed'],
        'pybind11.include.pybind11': ['*.h'],
        'pybind11.include.pybind11.detail': ['*.h'],
        'pybind11.include.pybind11.eigen': ['*.h'],
        'pybind11.include.pybind11.stl': ['*.h'],
        'pybind11.include.pybind11.conduit': ['*.h'],
        'pybind11.share.cmake.pybind11': ['*.cmake'],
        'pybind11.share.pkgconfig': ['*.pc'],
    },
    entry_points={
        'console_scripts': ['pybind11-config = pybind11.__main__:main'],
    },
)
" #$version))))))))
    (native-inputs
     (list cmake-minimal
           python-setuptools
           python-wheel
           python-wrapper))
    (home-page "https://github.com/pybind/pybind11/")
    (synopsis "Seamless operability between C++11 and Python")
    (description
     "@code{pybind11} is a lightweight header-only library that exposes C++
types in Python and vice versa, mainly to create Python bindings of existing
C++ code.  Its goals and syntax are similar to the @code{Boost.Python}
library: to minimize boilerplate code in traditional extension modules by
inferring type information using compile-time introspection.")
    (license license:bsd-3)))

(packages->manifest (list pybind11-3))
