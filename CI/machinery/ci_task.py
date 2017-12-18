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

#Some builds take really a long time
@timeout(10000)
def call(*args, **kwargs):
    return check_call(*args, **kwargs)


class CiTask():

    def __init__(self,
                 mode='Continuous',
                 build_configuration='Release',
                 distrib=None,
                 ci_config=None,
                 fast=True,
                 pkgs=None,
                 srcs=None,
                 targets=None,
                 cmake_cmd=None,
                 directories=[]):
        """Create a task, see examples in tasks.py.
        """

        self._fast = fast
        self._distrib = distrib
        self._mode = mode
        self._build_configuration = build_configuration
        self._ci_config = ci_config
        self._pkgs = pkgs
        self._srcs = srcs
        self._targets = targets
        self._cmake_cmd = cmake_cmd
        self._directories = directories

    def build_dir(self, src):
        """Return a name
        depending on src, distrib and task name
        """
        if isinstance(self._ci_config, str):
            ci_config_name = self._ci_config
        else:
            ci_config_name = '-'.join(self._ci_config)
        return ('build-' + src.replace('.', '-') + self._distrib.replace(':', '-') + '-' +
                ci_config_name).replace('/', '-')

    def templates(self):
        # remove build-base, gnu-c++, gfortran, it is redundant
        return ','.join(self._pkgs)

    def copy(self):
        """Copy constructor :
         duplicate a task and possibly extend configuration of the result.
        """
        def init(mode=self._mode,
                 build_configuration=self._build_configuration,
                 distrib=self._distrib,
                 ci_config=self._ci_config, fast=self._fast, pkgs=self._pkgs,
                 srcs=self._srcs, targets=self._targets,
                 cmake_cmd=self._cmake_cmd,
                 add_directories=None,
                 add_pkgs=None, remove_pkgs=None, add_srcs=None,
                 remove_srcs=None, add_targets=None, remove_targets=None):

            # WARNING: remember that default arg are mutable in python
            # http://docs.python-guide.org/en/latest/writing/gotchas/

            if add_pkgs is not None:
                pkgs = self._pkgs + add_pkgs

            if remove_pkgs is not None:
                pkgs = list(filter(lambda p: p not in remove_pkgs, pkgs))

            if add_srcs is not None:
                srcs = self._srcs + add_srcs

            if remove_srcs is not None:
                srcs = list(filter(lambda p: p not in remove_srcs, srcs))

            if add_targets is not None:
                targets = self._targets + add_targets

            if remove_targets is not None:
                targets = list(
                    filter(lambda p: p not in remove_targets, targets))

            if add_directories is not None:
                directories = self._directories + add_directories
            else:
                directories = self._directories

            return CiTask(mode, build_configuration, distrib, ci_config, fast,
                          pkgs, srcs, targets, cmake_cmd, directories)
        return init

    def run(self, root_dir, targets_override=None, dry_run=False):

        return_code = 0

        for src in self._srcs:
            # --- Path to CMakeLists.txt ---
            full_src = os.path.join(root_dir, src)
            if targets_override is not None:
                self._targets[src] = targets_override

            # --- Create build dir for src config ---
            bdir = self.build_dir(src)

            if not dry_run:
                if os.path.exists(bdir):
                    shutil.rmtree(bdir, ignore_errors=True)
                os.makedirs(bdir)


            #
            # not generic, to be moved somewhere else
            #
            redundants = [
                'build-base', 'gfortran', 'gnu-c++', 'lpsolve', 'wget', 'xz',
                'asan', 'cppunit_clang', 'python-env', 'profiling',
                'python3-env', 'path', 'h5py3']
            templ_list = [p.replace('+', 'x')
                          for p in self._pkgs if p not in redundants]

            # hm python is so lovely
            if isinstance(self._ci_config, str):
                ci_config_args = self._ci_config
            else:
                ci_config_args = ','.join(self._ci_config)

            # --- List of arguments for cmake command ---
            cmake_args = ['-DMODE={0}'.format(self._mode),
                          '-DCI_CONFIG={0}'.format(ci_config_args),
                          '-DWITH_DOCKER=1',
                          '-DBUILD_CONFIGURATION={0}'.format(
                              self._build_configuration),
                          '-DDOCKER_DISTRIB={0}'.format(self._distrib),
                          '-DDOCKER_TEMPLATES={0}'.format(self.templates()),
                          '-DDOCKER_TEMPLATE={0}'.format('-'.join(templ_list)),
                          '-DDOCKER_PROJECT_SOURCE_DIR={0}'.format(full_src)]

            if self._directories is not None:
                cmake_args.append('-DDOCKER_SHARED_DIRECTORIES={0}'.format(
                    ';'.join(self._directories)))

            # for examples ..
            if not os.path.samefile(root_dir, full_src):
                cmake_args.append('-DDOCKER_SHARED_DIRECTORIES={:}'.format(
                    root_dir))

            #
            # not generic, to be moved somewhere else
            #
            # probably broken
            if self._cmake_cmd:
                src_absolute_dir = os.path.normpath(
                    os.path.abspath(__file__) + '../../../..')
                cmake_args += [
                    '-DDOCKER_CMAKE_WRAPPER={:}/{:}'.format(src_absolute_dir,
                                                            self._cmake_cmd)]

            try:
                if os.path.exists(os.path.join(full_src, 'CI')):
                    full_cmd = ['cmake'] + cmake_args + [os.path.join(full_src,
                                                                      'CI')]
                else:
                    full_cmd = ['cmake'] + cmake_args + [full_src]
                if not dry_run:
                    print("cmake command is: {:}".format(' '.join(full_cmd)))
                    call(full_cmd, cwd=bdir)
                    for target in self._targets[src]:
                        call(['make'] + ['-ki'] + [target], cwd=bdir)
                else:
                    msg = 'Would call: \n  - {:}'.format(' '.join(full_cmd))
                    msg += '\n  - make - ki target, \n for target in '
                    msg += '{:}'.format(' '.join(self._targets[src]))
                    msg += '\n both from path ' + bdir
                    print msg

            except CalledProcessError as error:
                return_code = 1
                print(error)

        return return_code

    def clean(self):

        for src in self._srcs:

            bdir = self.build_dir(src)

            try:
                call(['make'] + ['docker-clean-usr-local'], cwd=bdir)

                if not self._fast:

                    call(['make'] + ['docker-clean'],
                         cwd=bdir)

            except CalledProcessError as error:
                    print(error)

    def display(self):
        for attr, value in self.__dict__.items():
            print('\t{:} = {:}'.format(attr, value))
