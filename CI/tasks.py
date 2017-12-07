""" List of all possible 'tasks', i.e. ci configurations.

WARNING : this is siconos specific!!

A task, see :class:`machinery.ci_tasks.CiTask` must be defined with at least:
* a name (ci_config)
* a distribution (name:version)
* a list of dependencies (pkgs)

"""
from machinery.ci_task import CiTask
import os

# PLEASE KEEP CONFIG AS WHAT THEY MEAN.
# DO NOT ADD PACKAGES IF THEY ARE NOT NECESSARY.

#
# 1. where the packages configurations are defined
# Used in driver.py.
database = os.path.join('config', 'siconos.yml')

#
# 2. the default task
#
default = CiTask(
    ci_config='default',
    distrib='ubuntu:16.04',
    pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++', 'atlas-lapack',
          'python-env'],
    srcs=['.'],
    targets={'.': ['docker-build', 'docker-ctest']})

minimal = CiTask(
    ci_config='minimal',
    distrib='ubuntu:16.10',
    pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++',
          'atlas-lapack', 'python-minimal'],
    srcs=['.'],
    targets={'.': ['docker-build', 'docker-ctest']})

minimal_with_python = CiTask(
    ci_config='minimal_with_python',
    distrib='ubuntu:16.10',
    pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++',
          'atlas-lapack', 'python-env'],
    srcs=['.'],
    targets={'.': ['docker-build', 'docker-ctest']})

#
# 3. all the tasks
#

siconos_default = default

siconos_default_nix = default.copy()(
    ci_config='nix',
    distrib='nixos/nix:latest',
    targets={'.': ['docker-build', 'docker-ctest']})

siconos_with_lpsolve = siconos_default.copy()(
    add_pkgs=['lpsolve'])

siconos_debian_latest = siconos_default.copy()(
    ci_config='with_bullet',
    add_pkgs=['bullet', 'h5py'],  # for mechanics.io
    distrib='debian:latest')

siconos_ubuntu_15_04 = siconos_default.copy()(
    distrib='ubuntu:15.04')

siconos_ubuntu_14_04 = siconos_default.copy()(
    distrib='ubuntu:14.04')

siconos_ubuntu_16_10 = siconos_default.copy()(
    distrib='ubuntu:16.10')

siconos_ubuntu_17_04 = siconos_default.copy()(
    distrib='ubuntu:17.04')

siconos_ubuntu_15_10 = siconos_default.copy()(
    distrib='ubuntu:15.10')

siconos_cxx_11_ubuntu_17_04 = siconos_default.copy()(
    distrib='ubuntu:17.04',
    ci_config='with_cxx11')


import os
from os.path import expanduser
home = expanduser("~")

siconos_documentation = siconos_default.copy()(
    distrib='ubuntu:16.10',
    ci_config='with_documentation',
    add_pkgs=['documentation'],
    add_directories=[os.path.join(home, '.ssh:/root/.ssh')],
    targets={'.': ['docker-build', 'docker-cmake', 'docker-make',
                   'docker-make-install',
                   'docker-make-doc', 'docker-make-upload']})

siconos_ubuntu_15_10_with_mechanisms = siconos_default.copy()(
    ci_config='with_mechanisms_conda_version',
    add_pkgs=['pythonocc-conda', 'wget', 'bash', 'bzip2',
              'pythonocc-conda-dep'],
    cmake_cmd='Build/ci-scripts/conda.sh',
    distrib='debian:stretch')

siconos_debian_mechanisms = siconos_default.copy()(
    ci_config='with_mechanisms',
    add_pkgs=['wget', 'bash', 'bullet', 'h5py', 'oce-pythonocc-deps',
              'oce-pythonocc'],
    distrib='debian:latest')


siconos_ubuntu_latest_mechanisms = siconos_default.copy()(
    ci_config='with_mechanisms',
    add_pkgs=['wget', 'bash', 'bullet', 'h5py', 'oce-pythonocc-deps',
              'oce-pythonocc'],
    distrib='ubuntu:latest')

siconos_numerics_only = siconos_ubuntu_17_04.copy()(
    ci_config='no_cxx',
    remove_pkgs=['gnu-c++'])

siconos_profiling = siconos_ubuntu_17_04.copy()(
    build_configuration='Profiling',
    add_pkgs=['profiling'])

# note fedora/atlas-lapack in siconos.yml -> cmake does not detect blas
siconos_fedora_latest = siconos_default.copy()(
    distrib='fedora:latest',
    ci_config=('with_umfpack',),
    remove_pkgs=['atlas-lapack', 'python-env'],
    add_pkgs=['openblas-lapacke', 'python3-env', 'umfpack'])

siconos_openblas_lapacke = siconos_default.copy()(
    remove_pkgs=['atlas-lapack'],
    add_pkgs=['openblas-lapacke'])

siconos_clang = siconos_ubuntu_17_04.copy()(
    ci_config=('with_bullet', 'with_py3'),
    remove_pkgs=['python-env'],
    add_pkgs=['clang-3.9', 'bullet', 'cppunit_clang-3.9', 'wget', 'xz', 'python3-env', 'path', 'h5py3'])  # h5py-3 for mechanics.io

siconos_clang_asan = siconos_clang.copy()(
    ci_config=('with_asan_clang', 'with_mumps', 'with_hdf5', 'with_serialization', 'with_py3'),
    add_pkgs=['mumps', 'hdf5', 'serialization'],
    build_configuration='Debug',)

# <clang-3.7.1 does not support linux 4.2
# This will likely hurt you
siconos_clang_msan = siconos_default.copy()(
    distrib='debian:jessie',
    ci_config='with_msan',
    build_configuration='Debug',
    add_pkgs=['clang-3.8', 'libcxx_msan', 'wget', 'xz', 'path'])

siconos_clang_cfi = siconos_default.copy()(
    distrib='debian:jessie',
    ci_config='with_cfi',
    build_configuration='Debug',
    add_pkgs=['mumps', 'hdf5', 'cfi'])

siconos_gcc_asan = siconos_fedora_latest.copy()(
    ci_config=('with_asan', 'with_mumps', 'with_hdf5', 'with_serialization'),
    cmake_cmd='Build/ci-scripts/fedora-mpi.sh',
    add_pkgs=['mumps', 'hdf5', 'asan', 'serialization', 'path', 'wget'],   # wget for path
    build_configuration='Debug')

siconos_gcc_asan_latest = siconos_fedora_latest.copy()(
    ci_config=('with_asan', 'with_mumps', 'with_hdf5', 'with_serialization'),
    distrib='fedora:rawhide',
    cmake_cmd='Build/ci-scripts/fedora-mpi.sh',
    add_pkgs=['mumps', 'hdf5', 'asan', 'serialization', 'path', 'wget'],   # wget for path
    build_configuration='Debug',
    fast=False)

# There is a bug in boost 1.58 distributed with Xenial (Ubuntu LTS 16.04).
# As long as it is not patched, we have to build on a newer ubuntu
siconos_serialization = siconos_ubuntu_17_04.copy()(
    ci_config='with_serialization',
    add_pkgs=['serialization'])

siconos_with_mumps = siconos_default.copy()(
    ci_config='with_mumps',
    add_pkgs=['mumps'])

siconos_with_umfpack = siconos_default.copy()(
    ci_config='with_umfpack',
    add_pkgs=['umfpack'])


# --- Config to run siconos examples ---

# Note FP/MB : this should be the only task(s) that run examples !!!
#        We may add later some 'examples-with-bullet' or 'examples-with-mechanisms' ... tasks later.
#        All tasks 'examples' should:
#         - conf, make and make install of siconos components (no tests!)
#         - conf, make and make test of siconos examples

# Case1 : siconos 'basics' components, numerics, kernel, control and related examples
siconos_light_examples = minimal_with_python.copy()(
    ci_config='examples_light',
    targets={'.': ['docker-build', 'docker-cmake', 'docker-make',
                   'docker-make-install', 'docker-make-clean'],
             'examples': ['docker-build', 'docker-ctest', 'docker-make-clean']},
    add_srcs=['examples'])

# Case2 : siconos with mechanics components and bullet + related examples
siconos_all_examples = minimal_with_python.copy()(
    ci_config='examples_all',
    add_pkgs=['bullet', 'h5py'],
    targets={'.': ['docker-build', 'docker-cmake', 'docker-make',
                   'docker-make-install', 'docker-make-clean'],
             'examples': ['docker-build', 'docker-ctest', 'docker-make-clean']},
    add_srcs=['examples'])

siconos_test_deb = CiTask(
    ci_config='examples',
    distrib='ubuntu:16.04',
    pkgs=['siconos'],
    srcs=['examples'])

siconos_test_rpm = CiTask(
    ci_config='examples',
    distrib='fedora:latest',
    pkgs=['siconos'],
    srcs=['examples'])

siconos_dev_mode_strict = siconos_default.copy()(
    ci_config='with_dev_mode_strict')

siconos_frama_c = siconos_default.copy()(
    ci_config='with_frama_c',
    add_pkgs=['opam', 'frama-c', 'libgtksourceview2.0-dev',
              'libgnomecanvas2-dev', 'aspcud', 'm4',
              'unzip', 'coq', 'ocaml', 'z3'])

#
# 4. dispatch based on hostname and distrib type (to min. disk requirement)
#
known_tasks = {'siconos---vm0':
               (siconos_fedora_latest,
                siconos_gcc_asan,
                siconos_gcc_asan_latest,
                siconos_debian_mechanisms,
                siconos_ubuntu_15_10),

               'siconos---vm1':
               (minimal,
                minimal_with_python,
                siconos_with_lpsolve,
                siconos_documentation,
                siconos_dev_mode_strict,
                siconos_clang,
                siconos_clang_asan),

               'siconos---vm2':
               (siconos_clang_msan,
                siconos_default_nix,
                siconos_ubuntu_15_10_with_mechanisms,
                siconos_ubuntu_15_04,
                siconos_ubuntu_14_04),

               'siconos---vm3':
               (siconos_debian_latest,
                siconos_openblas_lapacke,
                siconos_with_mumps,
                siconos_with_umfpack,
                siconos_light_examples,
                siconos_all_examples),

               'siconos---vm4':
               (siconos_profiling,
                siconos_ubuntu_17_04,
                siconos_numerics_only,
                siconos_cxx_11_ubuntu_17_04,
                siconos_serialization)}
