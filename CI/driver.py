#!/usr/bin/env python

from socket import gethostname

from tasks import siconos_default, known_tasks
from subprocess import check_output

hostname = gethostname().split('.')[0]

if hostname in known_tasks:
    tasks = known_tasks[hostname]
else:
    tasks = [siconos_default]

return_code = 0
for task in tasks:

    return_code += task.run()

for task in tasks:

    task.clean()

# clean everything (-> maybe once a week?)
def mklist(sstr):
    return filter(lambda s: s!='', sstr.strip().split('\n'))

running_containers=mklist(check_output(['docker', 'ps', '-q']))
if len(running_containers)>0:
    check_call(['docker', 'kill'] + running_containers)

containers=mklist(check_output(['docker', 'ps', '-a', '-q']))
if len(containers)>0:
    check_call(['docker', 'rm'] + containers)

images=mklist(check_output(['docker', 'images', '-q']))[1:]
if len(images)>0:
    check_call(['docker', 'rmi'] + images)

exit(return_code)
