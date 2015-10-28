from machinery.ci_task import CiTask

siconos_default = CiTask(
    ci_config='default',
    distrib='ubuntu:14.04',
    pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++', 'atlas-lapack'])

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

# dispatch based on hostname
known_tasks = {'fedora18-x86-64':
               [siconos_fedora_latest,
                siconos_openblas_lapacke,
                siconos_with_mumps],

               'siconos---fedora20':
               [siconos_clang,
                siconos_serialization],

               'siconos---ubuntu-12-04-amd64': 
               [siconos_default_profiling],

               'siconos---vm0':
               [siconos_ubuntu_15_10,
                siconos_ubuntu_15_04, 
                siconos_ubuntu_14_10]}
