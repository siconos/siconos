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

siconos_clang = siconos_ubuntu_15_10.copy()(
    add_pkgs=['clang'])

siconos_clang_asan = siconos_clang.copy()(
    ci_config='with_asan',
    add_pkgs=['mumps', 'hdf5'])

siconos_clang_msan = siconos_clang.copy()(
    ci_config='with_msan')

siconos_gcc_asan = siconos_openblas_lapacke.copy()(
    ci_config='with_asan',
    add_pkgs=['mumps', 'hdf5'])

siconos_serialization = siconos_default.copy()(
    ci_config='with_serialization',
    add_pkgs=['serialization'])

siconos_with_mumps = siconos_default.copy()(
    ci_config='with_mumps',
    add_pkgs=['mumps'])


siconos_default_examples = siconos_default.copy()(
    ci_config='examples',
    srcs=['.', 'examples'],
    targets={'.': ['docker-build', 'docker-cmake', 'docker-make',
                   'docker-make-install'],
             'examples': ['docker-build', 'docker-ctest']})


# dispatch based on hostname
known_tasks = {'siconos---vm0':
               [siconos_fedora_latest,
                siconos_openblas_lapacke,
                siconos_with_mumps,
                siconos_gcc_asan,
                siconos_debian_latest],

               'siconos---vm1':
               [siconos_default_examples,
                siconos_clang,
                siconos_serialization,
                siconos_clang_asan],
               #  siconos_clang_msan], not ready yet

               'siconos---vm2':
               [siconos_ubuntu_15_10,
                siconos_ubuntu_15_04,
                siconos_ubuntu_14_10,
                siconos_default_profiling]}
