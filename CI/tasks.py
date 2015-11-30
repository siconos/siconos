import os
from machinery.ci_task import CiTask

siconos_default = CiTask(
    ci_config='default',
    distrib='ubuntu:14.04',
    pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++', 'atlas-lapack'],
    srcs=['.'],
    targets={'.': ['docker-build', 'docker-ctest']})

siconos_default_profiling = siconos_default.copy()(
    build_configuration='Profiling')

siconos_debian_latest = siconos_default.copy()(
    distrib='debian:latest')

siconos_ubuntu_14_10 = siconos_default.copy()(
    distrib='ubuntu:14.10')

siconos_ubuntu_15_04 = siconos_default.copy()(
    distrib='ubuntu:15.04')

siconos_ubuntu_15_10 = siconos_default.copy()(
    distrib='ubuntu:15.10')

# note fedora/atlas-lapack in siconos.yml -> cmake does not detect blas
siconos_fedora_latest = siconos_default.copy()(
    distrib='fedora:latest',
    remove_pkgs=['atlas-lapack'],
    add_pkgs=['openblas-lapacke'])

siconos_openblas_lapacke = siconos_default.copy()(
    remove_pkgs=['atlas_lapack'],
    add_pkgs=['openblas-lapacke'])

siconos_clang = siconos_default.copy()(
    add_pkgs=['clang'],
    remove_pkgs=['gcc', 'gnu-c++'])

siconos_serialization = siconos_default.copy()(
    ci_config='with_serialization',
    add_pkgs=['serialization'])

siconos_with_mumps = siconos_default.copy()(
    ci_config='with_mumps',
    add_pkgs=['mumps'])

siconos_default_examples = siconos_default.copy()(
    ci_config='examples',
    add_srcs=['examples'],
    targets={'.': ['docker-build', 'docker-cmake', 'docker-make', 'docker-make-install'],
             'examples': ['docker-build', 'docker-ctest']})


# dispatch based on hostname
known_tasks = {'siconos---vm0':
               [siconos_fedora_latest,
                siconos_openblas_lapacke,
                siconos_with_mumps],

               'siconos---vm1':
               [siconos_clang,
                siconos_serialization],

               'siconos---ubuntu-12-04-amd64': 
               [siconos_default_examples,
                siconos_default_profiling],

               'siconos---vm2':
               [siconos_ubuntu_15_10,
                siconos_ubuntu_15_04, 
                siconos_ubuntu_14_10]}
