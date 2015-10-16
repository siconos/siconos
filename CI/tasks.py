from machinery.ci_task import CiTask

siconos_default = CiTask()

siconos_openblas_lapacke = CiTask(pkgs=['build-base','gcc', 'gfortran', 'g++',
                                        'openblas-lapacke'])

siconos_clang = CiTask(pkgs=['build-base','clang', 'gfortran',
                             'atlas-lapack'])

siconos_serialization = CiTask(ci_config='with_serialization',
                               pkgs = siconos_default._pkgs + ['serialization'])
siconos_debian_latest = CiTask(distrib='debian:latest')
siconos_fedora_latest = CiTask(distrib='fedora:latest')
siconos_with_mumps = CiTask(pkgs = siconos_default._pkgs + ['mumps'],
                            ci_config='with_mumps')

# dispatch based on hostname
known_tasks = {'fedora18-x86-64': [siconos_with_mumps,
                                   siconos_openblas_lapacke,
                                   siconos_clang,
                                   siconos_serialization]}
