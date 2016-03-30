from machinery.ci_task import CiTask

siconos_default = CiTask(
    ci_config='default',
    distrib='ubuntu:14.04',
    pkgs=['build-base', 'gcc', 'gfortran', 'gnu-c++', 'atlas-lapack', 'lpsolve', 'python-env'],
    srcs=['.'],
    targets={'.': ['docker-build', 'docker-ctest']})

siconos_test_deb = CiTask(
    ci_config='examples',
    distrib='ubuntu:15.10',
    pkgs=['siconos'],
    srcs=['examples'],
    targets={'examples': ['docker-build', 'docker-ctest']})

siconos_test_rpm = CiTask(
    ci_config='examples',
    distrib='fedora:latest',
    pkgs=['siconos'],
    srcs=['examples'],
    targets={'examples': ['docker-build', 'docker-ctest']})

siconos_debian_latest = siconos_default.copy()(
    ci_config='with_bullet',
    add_pkgs=['bullet', 'h5py'],  # for mechanics.io
    distrib='debian:latest')

siconos_ubuntu_15_04 = siconos_default.copy()(
    distrib='ubuntu:15.04')

siconos_ubuntu_15_10 = siconos_default.copy()(
    ci_config='with_umfpack',
    add_pkgs=['umfpack'],
    distrib='ubuntu:15.10')

siconos_ubuntu_15_10_with_mechanisms = siconos_default.copy()(
    ci_config='with_mechanisms',
    add_pkgs=['pythonocc-conda', 'wget', 'bash', 'bzip2', 'pythonocc-conda-dep'],
    cmake_cmd='Build/ci-scripts/conda.sh',
    distrib='debian:stretch')

siconos_numerics_only = siconos_ubuntu_15_10.copy()(
    ci_config='no_cxx',
    remove_pkgs=['gnu-c++'])

siconos_profiling = siconos_ubuntu_15_10.copy()(
    build_configuration='Profiling',
    add_pkgs=['profiling'])

# note fedora/atlas-lapack in siconos.yml -> cmake does not detect blas
siconos_fedora_latest = siconos_default.copy()(
    distrib='fedora:latest',
    ci_config=('with_umfpack',),
    remove_pkgs=['atlas-lapack', 'python-env'],
    add_pkgs=['openblas-lapacke', 'python3-env', 'umfpack'])

siconos_openblas_lapacke = siconos_default.copy()(
    ci_config='with_umfpack',
    remove_pkgs=['atlas-lapack'],
    add_pkgs=['openblas-lapacke', 'umfpack', 'path', 'wget'],  # wget for path
    with_examples=True)

siconos_clang = siconos_ubuntu_15_10.copy()(
    ci_config=('with_bullet', 'with_py3'),
    with_examples=False,
    remove_pkgs=['python-env'],
    add_pkgs=['clang', 'bullet', 'cppunit_clang', 'wget', 'xz', 'python3-env', 'path', 'h5py3'])  # h5py-3 for mechanics.io

siconos_clang_asan = siconos_clang.copy()(
    ci_config=('with_asan_clang', 'with_mumps', 'with_hdf5', 'with_serialization', 'with_py3'),
    add_pkgs=['mumps', 'hdf5', 'serialization'],
    build_configuration='Debug',
    with_examples=False)

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
    add_pkgs=['clang-3.8', 'mumps', 'hdf5', 'cfi'])

siconos_gcc_asan = siconos_fedora_latest.copy()(
    ci_config=('with_asan', 'with_mumps', 'with_hdf5', 'with_serialization'),
    cmake_cmd='Build/ci-scripts/fedora-mpi.sh',
    add_pkgs=['mumps', 'hdf5', 'asan', 'serialization', 'path', 'wget'],   # wget for path
    build_configuration='Debug',
    with_examples=False)

siconos_gcc_asan_latest = siconos_fedora_latest.copy()(
    ci_config=('with_asan', 'with_mumps', 'with_hdf5', 'with_serialization'),
    distrib='fedora:rawhide',
    cmake_cmd='Build/ci-scripts/fedora-mpi.sh',
    add_pkgs=['mumps', 'hdf5', 'asan', 'serialization', 'path', 'wget'],   # wget for path
    build_configuration='Debug',
    fast=False)

siconos_serialization = siconos_default.copy()(
    ci_config='with_serialization',
    add_pkgs=['serialization'])

siconos_with_mumps = siconos_default.copy()(
    ci_config='with_mumps',
    add_pkgs=['mumps'])


siconos_default_examples = siconos_default.copy()(
    ci_config='examples',
    with_examples=True)


# dispatch based on hostname and distrib type (to min. disk requirement)

known_tasks = {'siconos---vm0':
               (siconos_fedora_latest,
                siconos_gcc_asan,
                siconos_gcc_asan_latest,
                siconos_ubuntu_15_10_with_mechanisms),

               'siconos---vm1':
               (siconos_numerics_only,
                siconos_clang,
                siconos_clang_asan,
                siconos_clang_msan),

               'siconos---vm2':
               (siconos_ubuntu_15_10,
                siconos_ubuntu_15_04,
                siconos_profiling),

               'siconos---vm3':
               (siconos_debian_latest,
                siconos_openblas_lapacke,
                siconos_serialization,
                siconos_with_mumps)}
