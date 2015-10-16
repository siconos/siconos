#!/usr/bin/env python

import os
from subprocess import check_call, CalledProcessError
from socket import gethostname

from tasks import siconos_default, known_tasks

hostname = gethostname().split('.')[0]

if hostname in known_tasks:
    tasks = known_tasks[hostname]
else:
    tasks = [siconos_default]

return_code = 0
for task in tasks:

    if not os.path.exists(task._ci_config):
        os.makedirs(task._ci_config)

    cmake_args = ['-DMODE={0}'.format(task._mode),
                  '-DCI_CONFIG={0}'.format(task._ci_config),
                  '-DWITH_DOCKER=1',
                  '-DDOCKER_DISTRIB={0}'.format(task._distrib),
                  '-DDOCKER_TEMPLATES={0}'.format(task.templates())]

    try:
        check_call(['cmake'] + cmake_args + [os.path.join('..', '..')],
                   cwd=task._ci_config)
        check_call(['make'] + ['docker-build'],
                   cwd=task._ci_config)
        check_call(['make'] + ['docker-ctest'],
                   cwd=task._ci_config)

        if not task._fast:
            check_call(['make'] + ['docker-clean'],
                       cwd=task._ci_config)

    except CalledProcessError as error:
        return_code = 1
        print error

exit(return_code)
