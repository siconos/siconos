from machinery.ci_task import CiTask

siconos_default = CiTask(
    ci_config='default',
    distrib='ubuntu:14.04',
    pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++', 'atlas-lapack'])

siconos_debian_latest = siconos_default.copy()(
    distrib='debian:latest')

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
known_tasks = {'fedora18-x86-64': [siconos_fedora_latest,
                                   siconos_openblas_lapacke,
                                   siconos_clang,
                                   siconos_with_mumps,
                                   siconos_serialization]}
