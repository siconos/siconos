from machinery.ci_task import CiTask

siconos_default = CiTask()

siconos_serialization = CiTask(ci_config='serialization')
siconos_debian_latest = CiTask(distrib='debian:latest')
siconos_fedora_latest = CiTask(distrib='fedora:latest')

# dispatch based on hostname
known_tasks = {'fedora18-x86_64': [siconos_default,
                                   siconos_serialization]}
