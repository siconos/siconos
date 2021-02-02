""" List of all possible 'tasks', i.e. ci configurations.

WARNING : this is siconos specific!!

A task, see :class:`machinery.ci_tasks.CiTask` must be defined with at least:
* a name (ci_config)
* a distribution (name:version)
* a list of dependencies (pkgs)

"""
from machinery.ci_task import SiconosCiTask
from os.path import expanduser

# PLEASE KEEP CONFIGS AS WHAT THEY MEAN.
# DO NOT ADD PACKAGES IF THEY ARE NOT NECESSARY.

empty = SiconosCiTask()

# 1. Define a default task, that will be used as basis for all other tasks.
# - build a docker container
# - distrib : docker image name
# - ci_config : siconos config(s) name(s)
# - pkgs : list of deps to be installed
# - srcs : path to main CMakeLists.txt
# - targets : list of cmake targets to be executed.
# - fast : default = True, set false to clean properly docker containers before
#          run.
siconos_default = SiconosCiTask(
    docker=True,
    ci_config='default',
    distrib='ubuntu:18.04',
    pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++', 'openblas-lapacke',
          'python3-env'],
    srcs=['.'],
    targets={'.': ['docker-build', 'docker-ctest', 'docker-submit']})

default = siconos_default  # shortcut for .travis.yml

#
# 2. Define all the tasks
#

siconos_default_nix = siconos_default.copy()(
    ci_config='nix',
    distrib='nixos/nix:latest')

siconos_debian_latest = siconos_default.copy()(
    distrib='debian:latest')

siconos_ubuntu_18_04 = siconos_default.copy()(
    distrib='ubuntu:18.04')

siconos_ubuntu_20_04 = siconos_default.copy()(
    distrib='ubuntu:20.04')

siconos_fedora_latest = siconos_default.copy()(
    distrib='fedora:latest')

siconos_gazebo = siconos_default.copy()(
    distrib='nvidia/opengl:1.0-glvnd-devel-ubuntu16.04',
    ci_config=('with_bullet', 'with_py3'),
    add_pkgs=['bullet', 'gazebo'],
    targets={'.': ['docker-build', 'docker-cmake', 'docker-make',
                   'docker-make-install', 'docker-cmd']})

siconos_with_lpsolve = siconos_default.copy()(
    add_pkgs=['lpsolve'])

home = expanduser("~")

siconos_debian_mechanisms = siconos_default.copy()(
    ci_config='with_mechanisms',
    add_pkgs=['wget', 'bash', 'h5py', 'oce-pythonocc-deps',
              'oce-pythonocc'],
    distrib='debian:latest')


siconos_ubuntu_latest_mechanisms = siconos_default.copy()(
    ci_config='with_mechanisms',
    add_pkgs=['wget', 'bash', 'h5py', 'oce-pythonocc-deps',
              'oce-pythonocc'],
    distrib='ubuntu:latest')

siconos_numerics_only = siconos_ubuntu_18_04.copy()(
    ci_config='no_cxx',
    remove_pkgs=['gnu-c++'])

siconos_profiling = siconos_ubuntu_18_04.copy()(
    build_configuration='Profiling',
    add_pkgs=['profiling'])

siconos_fedora_latest_with_umfpack = siconos_default.copy()(
    distrib='fedora:latest',
    ci_config=('with_umfpack',),
    remove_pkgs=['python-env'],
    add_pkgs=['python3-env', 'umfpack'])

siconos_clang = siconos_ubuntu_18_04.copy()(
    ci_config=('with_bullet', 'with_py3'),
    remove_pkgs=['python-env'],
    add_pkgs=['clang-3.9', 'bullet', 'cppunit_clang-3.9', 'wget', 'xz',
              'python3-env', 'path', 'h5py3'])  # h5py-3 for mechanics.io

siconos_clang_asan = siconos_clang.copy()(
    ci_config=('with_asan_clang', 'with_mumps', 'with_hdf5',
               'with_serialization', 'with_py3'),
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
    #    cmake_cmd='Build/ci-scripts/fedora-mpi.sh',
    add_pkgs=['mumps', 'hdf5', 'asan', 'serialization', 'path', 'wget'],
    # wget for path
    build_configuration='Debug')

siconos_gcc_asan_latest = siconos_fedora_latest.copy()(
    ci_config=('with_asan', 'with_mumps', 'with_hdf5', 'with_serialization'),
    distrib='fedora:rawhide',
    #    cmake_cmd='Build/ci-scripts/fedora-mpi.sh',
    add_pkgs=['mumps', 'hdf5', 'asan', 'serialization', 'path', 'wget'],
    # wget for path
    build_configuration='Debug',
    fast=False)

siconos_serialization = siconos_ubuntu_18_04.copy()(
    ci_config='with_serialization',
    add_pkgs=['serialization'])

siconos_with_mumps = siconos_default.copy()(
    ci_config='with_mumps',
    add_pkgs=['mumps'])

siconos_with_umfpack = siconos_default.copy()(
    ci_config='with_umfpack',
    add_pkgs=['umfpack'])


siconos_dev_mode_strict = siconos_default.copy()(
    ci_config='with_dev_mode_strict')

siconos_frama_c = siconos_default.copy()(
    ci_config='with_frama_c',
    add_pkgs=['opam', 'frama-c', 'libgtksourceview2.0-dev',
              'libgnomecanvas2-dev', 'aspcud', 'm4',
              'unzip', 'coq', 'ocaml', 'z3'])

