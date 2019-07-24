#!/usr/bin/env python

import sys
import os
from socket import gethostname

from tasks import default, known_tasks
from getopt import gnu_getopt, GetoptError

import functools

def usage():
    print("""
    {0} [-v] [--packages-db=...] \
        [--task=<task1>,<task2>,...] \
        [--list-tasks] \
        [--root-dir=...] \
        [--targets=...] \
        [--run]\
        [--brutal-docker-clean\
        [--print='script'|'docker'|'vagrant']\
        [--dry-run]
    """.format(sys.argv[0]))


try:
    opts, args = gnu_getopt(sys.argv[1:], 'v',
                            ['run', 'task=', 'list-tasks', 'targets=',
                             'packages-db=',
                             'root-dir=',
                             'print=', 'brutal-docker-clean', 'dry-run'])

except GetoptError as err:
    sys.stderr.write(str(err))
    usage()
    exit(2)

task = None
verbose = False
run = False
return_code = 0
print_mode = False
dry_run = False
database = None

for o, a in opts:
    if o in ('--run',):
        run = True

    if o in ('-v',):
        verbose = True

    if o in ('--packages-db',):
        database = a
        assert database is not None
        
    if o in ('--root-dir',):
        root_dir = os.path.abspath(a)

    if o in ('--task',):
        import tasks
        
        def arg_check(acc, kv):
            k, v = kv.split('=')
            return_val = None
            
            if k[-1] == '+':
                k = 'add_{0}'.format(str(k[:-1]))

            elif k[-1] == '-':
                k = 'remove_{0}'.format(str(k[:-1]))

            if v.lower() in ['false', 'off']:
                return_val = False
            
            elif v.lower() in ['true', 'on']:
                return_val = True

            elif ',' in v:
                return_val = list(filter(lambda s: s != '', v.split(',')))
                
            else:
                return_val = v

            if k in acc:
                if type(acc[k]) == list:
                    l = acc[k]
                    acc[k] = {l[0]:l[1:]}
                acc[k][return_val[0]] = return_val[1:]
            else:
                acc[k] = return_val

            return acc
                
        task_arg = a.split(':')
        task_name = task_arg[0]
        task_parameters = functools.reduce(arg_check, task_arg[1:], dict())
        task = getattr(tasks, task_name).copy()(**task_parameters)

    if o in ('--list-tasks',):
        import tasks
        from machinery.ci_task import CiTask
        for s in dir(tasks):
            if isinstance(getattr(tasks, s), CiTask):
                print(s)
                if verbose:
                    t = getattr(tasks, s)
                    t.display()
        sys.exit()
    
    if o in ('--print',):
        print_mode = True
        output_mode_str = a

    if o in ('--dry-run',):
        dry_run = True

        
hostname = gethostname().split('.')[0]

if task is None:

    if hostname in known_tasks:
        tasks = known_tasks[hostname]
    else:
        tasks = [default]
else:
    tasks = [task]

    
if print_mode:

    from machinery.mksenv import print_commands, output_mode_spec

    for task in tasks:
        # For each task, print the 'command' file (dockerfile, script ...)
        # for the given configuration
        distrib, distrib_version = task._distrib.split(':')
        print_commands(specfilename=os.path.join(os.path.dirname(os.path.realpath(__file__)), database),
                       distrib=distrib,
                       distrib_version=distrib_version,
                       pkgs=task._pkgs,
                       output_mode=output_mode_spec[output_mode_str],
                       split=False)

if run:

    for task in tasks:

        return_code += task.run(root_dir,
                                dry_run=dry_run)

        print('return code {0}'.format(return_code))

    if not dry_run:
        for task in tasks:
            try:
                return_code += task.clean()
                print('return code {0}'.format(return_code))

            except Exception as e:
                return_code += 1
                sys.stderr.write(str(e))

    exit(return_code)
