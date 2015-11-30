import os
from subprocess import check_call, CalledProcessError

class CiTask():

    def __init__(self,
                 mode='Continuous',
                 build_configuration='Release',
                 distrib=None,
                 ci_config=None,
                 fast=False,
                 pkgs=None,
                 srcs=None,
                 targets=None):
        self._fast = fast
        self._distrib = distrib
        self._mode = mode
        self._build_configuration = build_configuration
        self._ci_config = ci_config
        self._pkgs = pkgs
        self._srcs = srcs
        self._targets = targets


    def templates(self):
        return ','.join(self._pkgs)

    def copy(self):
        def init(mode=self._mode,
                 build_configuration=self._build_configuration,
                 distrib=self._distrib,
                 ci_config=self._ci_config, fast=self._fast, pkgs=self._pkgs,
                 srcs=self._srcs, targets=self._targets,
                 add_pkgs=None, remove_pkgs=None, add_srcs=None,
                 remove_srcs=None, add_targets=None, remove_targets=None):

            if add_pkgs is not None:
                pkgs = self._pkgs + add_pkgs

            if remove_pkgs is not None:
                pkgs = filter(lambda p: p not in remove_pkgs, pkgs)

            if add_srcs is not None:
                srcs += self._srcs + add_srcs

            if remove_srcs is not None:
                srcs = filter(lambda p: p not in remove_srcs, srcs)

            if add_targets is not None:
                targets += self._targets + add_targets

            if remove_targets is not None:
                targets = filter(lambda p: p not in remove_targets, targets)

            return CiTask(mode, build_configuration, distrib, ci_config, fast,
                          pkgs, srcs, targets)
        return init

    def run(self):

        return_code = 0

        for src in self._srcs:

            bdir = src.replace('.', '_')
            bdir += '_' + self._ci_config

            if not os.path.exists(bdir):
                os.makedirs(bdir)

            cmake_args = ['-DMODE={0}'.format(self._mode),
                          '-DCI_CONFIG={0}'.format(self._ci_config),
                          '-DWITH_DOCKER=1',
                          '-DBUILD_CONFIGURATION={0}'.format(
                              self._build_configuration),
                          '-DDOCKER_DISTRIB={0}'.format(self._distrib),
                          '-DDOCKER_TEMPLATES={0}'.format(self.templates())]

            try:
                check_call(['cmake'] + cmake_args + [os.path.join('..','..',src)],
                           cwd=bdir)

                for target in self._targets[src]:

                    check_call(['make'] + [target], cwd=bdir)

            except CalledProcessError as error:
                return_code = 1
                print error

        return return_code


    def clean(self):

        for src in self._srcs:

            bdir = src.replace('.', '_')
            bdir += '_' + self._ci_config
            if not self._fast:

                try:
                    check_call(['make'] + ['docker-clean'],
                               cwd=bdir)

                except CalledProcessError as error:
                    print error
