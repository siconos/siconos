"""CI task management

Note : this should remain independant of siconos.

"""
import os
import shutil
from subprocess import check_call, CalledProcessError
import time
import multiprocessing


class TimeoutException(Exception):
    pass


class RunableProcessing(multiprocessing.Process):
    def __init__(self, func, *args, **kwargs):
        self.queue = multiprocessing.Queue(maxsize=1)
        args = (func,) + args
        multiprocessing.Process.__init__(self, target=self.run_func, args=args,
                                         kwargs=kwargs)

    def run_func(self, func, *args, **kwargs):
        try:
            result = func(*args, **kwargs)
            self.queue.put((True, result))
        except Exception as e:
            self.queue.put((False, e))

    def done(self):
        return self.queue.full()

    def result(self):
        return self.queue.get()


def timeout(seconds, force_kill=True):
    if seconds == 0:
        def wrapper(function):
            return function
        return wrapper
    else:
        def wrapper(function):
            def inner(*args, **kwargs):
                now = time.time()
                proc = RunableProcessing(function, *args, **kwargs)
                proc.start()
                proc.join(seconds)
                if proc.is_alive():
                    if force_kill:
                        proc.terminate()
                    runtime = int(time.time() - now)
                    raise TimeoutException(
                        'timed out after {0} seconds'.format(runtime))
                if not proc.done():
                    proc.terminate()
                assert proc.done()
                success, result = proc.result()
                if success:
                    # return time.time() - now, result
                    return result
                else:
                    # raise time.time() - now, result
                    raise result
            return inner
        return wrapper


# Some builds take really a long time
@timeout(10000)
def call(*args, **kwargs):
    try:

        return_code = check_call(*args, **kwargs)

    except CalledProcessError as error:

        print(error)
        return_code = 1

    return return_code


class CiTask():

    def __init__(self,
                 mode='Continuous',
                 build_configuration='Release',
                 docker=False,
                 distrib=None,
                 ci_config=None,
                 fast=True,
                 pkgs=None,
                 srcs=None,
                 targets=None,
                 cmake_cmd='cmake',
                 make_cmd='make',
                 cmake_args=[],
                 directories=[]):

        """Create a task, see examples in tasks.py.
        """

        self._docker = docker
        self._fast = fast
        self._distrib = distrib
        self._mode = mode
        self._build_configuration = build_configuration
        self._ci_config = ci_config
        self._pkgs = pkgs
        self._srcs = srcs
        self._targets = targets
        self._cmake_cmd = cmake_cmd
        self._make_cmd = make_cmd
        self._cmake_args = cmake_args
        self._directories = directories

    def template_maker(self):
        assert False
        return '-'.join(self._pkgs)

    def build_dir(self, src):
        """Return a name
        depending on src, distrib and task name
        """
        if self._distrib is None:
            assert self._docker is False
            import platform
            distrib = platform.platform()
        else:
            distrib = self._distrib

        if isinstance(self._ci_config, str):
            ci_config_name = self._ci_config
        else:
            ci_config_name = '-'.join(self._ci_config)
        return ('build-' + src.replace('.', '-') + distrib.replace(':', '-') + '-' +
                ci_config_name).replace('/', '-')

    def templates(self):
        return ','.join(self._pkgs)

    def copy(self):
        """
         duplicate a task and possibly extend configuration of the result.
        """
        def init(mode=self._mode,
                 docker=self._docker,
                 build_configuration=self._build_configuration,
                 distrib=self._distrib,
                 ci_config=self._ci_config, fast=self._fast, pkgs=self._pkgs,
                 srcs=self._srcs, targets=self._targets,
                 cmake_cmd=self._cmake_cmd,
                 make_cmd=self._make_cmd,
                 cmake_args=self._cmake_args,
                 add_directories=None,
                 add_pkgs=None, remove_pkgs=None, add_srcs=None,
                 remove_srcs=None, add_targets=None, remove_targets=None):


            # WARNING: remember that default arg are mutable in python
            # http://docs.python-guide.org/en/latest/writing/gotchas/

            new_targets = dict()

            new_distrib = None

            if type(distrib) == list:
                new_distrib = ':'.join(distrib)
            else:
                if distrib is not None:
                    assert type(distrib) == str
                    new_distrib = distrib

            if type(targets) == list:
                for src in self._targets.keys():
                    new_targets[src] = targets

            else:
                assert type(targets) == dict
                new_targets = targets

            if add_pkgs is not None:
                pkgs = self._pkgs + add_pkgs

            if remove_pkgs is not None:
                pkgs = list(filter(lambda p: p not in remove_pkgs, pkgs))

            if add_srcs is not None:
                srcs = self._srcs + add_srcs

            if remove_srcs is not None:
                srcs = list(filter(lambda p: p not in remove_srcs, srcs))

            if add_targets is not None:
                for src in self._targets.keys():
                    new_targets[src] += add_targets

            if remove_targets is not None:
                for src in self._targets.keys():
                    new_targets[src] = list(
                        filter(lambda p: p not in remove_targets, self._targets[src]))

            if add_directories is not None:
                directories = self._directories + add_directories
            else:
                directories = self._directories

            from copy import deepcopy

            new_task = deepcopy(self)
            
            new_task.__init__(mode=mode, build_configuration=build_configuration,
                              docker=docker,
                              distrib=new_distrib, ci_config=ci_config, fast=fast,
                              pkgs=pkgs, srcs=srcs, targets=new_targets, cmake_cmd=cmake_cmd,
                              cmake_args=cmake_args,
                              make_cmd=make_cmd,
                              directories=directories)
            return new_task
            
        return init

    def run(self, root_dir, dry_run=False):

        return_code = 0

        for src in self._srcs:
            # --- Path to CMakeLists.txt ---
            full_src = os.path.join(root_dir, src)

            # --- Create build dir for src config ---
            bdir = self.build_dir(src)

            if not dry_run:
                if os.path.exists(bdir):
                    shutil.rmtree(bdir, ignore_errors=True)
                os.makedirs(bdir)


            # hm python is so lovely
            if isinstance(self._ci_config, str):
                ci_config_args = self._ci_config
            else:
                ci_config_args = ','.join(self._ci_config)

            # --- List of arguments for cmake command ---
            cmake_args = self._cmake_args
            if self._docker:
                cmake_args += ['-DMODE={0}'.format(self._mode),
                               '-DCI_CONFIG={0}'.format(ci_config_args),
                               '-DWITH_DOCKER=1',
                               '-DBUILD_CONFIGURATION={0}'.format(
                                   self._build_configuration),
                               '-DDOCKER_DISTRIB={0}'.format(self._distrib),
                               '-DDOCKER_TEMPLATES={0}'.format(self.templates()),
                               '-DDOCKER_TEMPLATE={0}'.format(self.template_maker()),
                               '-DDOCKER_PROJECT_SOURCE_DIR={0}'.format(full_src)]


            if self._docker and self._directories is not None:
                cmake_args.append('-DDOCKER_SHARED_DIRECTORIES={0}'.format(
                    ';'.join(self._directories)))

            # for examples ..
            if self._docker and not os.path.samefile(root_dir, full_src):
                cmake_args.append('-DDOCKER_SHARED_DIRECTORIES={:}'.format(
                    root_dir))

            try:
                if os.path.exists(os.path.join(full_src, 'CI')):
                    if self._docker:
                        full_cmd = [self._cmake_cmd] + cmake_args + [os.path.join(full_src,
                                                                          'CI')]
                    else:
                        full_cmd = [self._cmake_cmd] + cmake_args + [full_src]
                else:
                    full_cmd = [self._cmake_cmd] + cmake_args + [full_src]
                if not dry_run:
                    print("cmake command is: {:}".format(' '.join(full_cmd)))
                    return_code += call(full_cmd, cwd=bdir)
                    for target in self._targets[src]:
                        return_code += call([self._make_cmd] + [target], cwd=bdir)
                else:
                    msg = 'Would call: \n  - {:}'.format(' '.join(full_cmd))
                    msg += '\n  - make target, \n for target in '
                    msg += '{:}'.format(' '.join(self._targets[src]))
                    msg += '\n both from path ' + bdir
                    print (msg)

            except Exception as error:
                return_code = 1
                print(error)

        return return_code

    def clean(self):

        return_code = 0

        for src in self._srcs:

            bdir = self.build_dir(src)

            try:
                if self._docker:
                    return_code += call([self._make_cmd] + ['docker-clean-usr-local'], cwd=bdir)

                    if not self._fast:

                        return_code += call([self._make_cmd] + ['docker-clean'],
                                            cwd=bdir)

            except Exception as error:
                    print(error)

        return return_code

    def display(self):
        for attr, value in self.__dict__.items():
            print('\t{:} = {:}'.format(attr, value))
