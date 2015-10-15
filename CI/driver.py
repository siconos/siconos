#!/usr/bin/env python

from subprocess import check_call
from socket import gethostname

from tasks import siconos_default, known_tasks

if gethostname() in known_tasks:
   tasks = known_tasks[hostname]
else:
   tasks = [siconos_default]

for task in tasks:

   cmake_args = ['-DMODE={0}'.format(task._mode),
                 '-DCI_CONFIG={0}'.format(task._ci_config),
                 '-DWITH_DOCKER=1',
                 '-DDOCKER_DISTRIB={0}'.format(task._distrib),
                 '-DDOCKER_TEMPLATES={0}'.format(task.templates())]

   check_call(['cmake'] + cmake_args + ['..'])
   check_call(['make'] + ['docker-build'])
   check_call(['make'] + ['docker-ctest'])

