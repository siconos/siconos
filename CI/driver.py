#!/usr/bin/env python

import sys
from socket import gethostname

from tasks import default, known_tasks, database
from subprocess import check_output, check_call

from getopt import gnu_getopt, GetoptError


def usage():
    print("""
{0} [-v] [--tasks=<task1>,<task2>,...] [--list-tasks] [--targets] [--run]\
    [--brutal-docker-clean]
    """.format(sys.argv[0]))


try:
    opts, args = gnu_getopt(sys.argv[1:], 'v',
                            ['run', 'tasks=', 'list-tasks', 'targets=',
                             'print=', 'brutal-docker-clean'])

except GetoptError as err:
    sys.stderr.write(str(err))
    usage()
    exit(2)


tasks = None
brutal_clean = False
targets_override = None
verbose = False
run = False
return_code = 0
print_mode = False

for o, a in opts:
    if o in ('--run',):
        run = True

    if o in ('-v',):
        verbose = True

    if o in ('--tasks',):
        import tasks
        tasks = [getattr(tasks, s) for s in a.split(',')]

    if o in ('--list-tasks',):
        import tasks
        from machinery.ci_task import CiTask
        for s in dir(tasks):
            if isinstance(getattr(tasks, s), CiTask):
                print(s)
                if verbose:
                    t = getattr(tasks, s)
                    t.display()

    if o in ('--print',):
        print_mode = True
        output_mode_str = a

    if o in ('--targets',):
        targets_override = a.split(',')

    if o in ('--brutal-docker-clean',):
        brutal_clean = True

hostname = gethostname().split('.')[0]

if tasks is None:

    if hostname in known_tasks:
        tasks = known_tasks[hostname]
    else:
        tasks = [default]

if print_mode:

    from machinery.mksenv import print_commands, OutputMode, output_mode_spec

    for task in tasks:
        distrib, distrib_version = task._distrib.split(':')
        print_commands(specfilename=database,
                       distrib=distrib,
                       distrib_version=distrib_version,
                       pkgs=task._pkgs,
                       output_mode=output_mode_spec[output_mode_str],
                       split=False)
if run:

    for task in tasks:

        return_code += task.run(targets_override=targets_override)

    for task in tasks:

        task.clean()

# docker specific, to be moved elsewhere
if brutal_clean:

    # clean everything (-> maybe once a week?)
    def mklist(sstr):
        return filter(lambda s: s != '', sstr.strip().split('\n'))

    running_containers = mklist(check_output(['docker', 'ps', '-q']))

    if len(running_containers) > 0:
        check_call(['docker', 'kill'] + running_containers)

    containers = mklist(check_output(['docker', 'ps', '-a', '-q']))
    if len(containers) > 0:
        check_call(['docker', 'rm'] + containers)

    images = mklist(check_output(['docker', 'images', '-q']))[1:]
    if len(images) > 0:
        check_call(['docker', 'rmi'] + images)

    exit(return_code)
