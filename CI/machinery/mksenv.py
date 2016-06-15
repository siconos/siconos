#!/usr/bin/env python

#
# Build shell commands from package specifications in a yml database.
#

import getopt
import yaml
import sys
import shlex


def usage():
    print("""
{0} [--pkg=<pkg>] [--pkgs=<pkg1,pkg2>,...] [--script] \
[--docker] [--vagrant] [--split=...] [--distrib=...] \
    /path/to/<example>.yml""".format(sys.argv[0]))


class OutputMode:
    Script, Docker, Vagrant = range(3)

output_mode_spec = dict()
output_mode_spec['script'] = OutputMode.Script
output_mode_spec['docker'] = OutputMode.Docker
output_mode_spec['vagrant'] = OutputMode.Vagrant

def wildcard(spec):
    if 'wildcard' in spec:
        return spec['wildcard']
    else:
        return 'any'


def is_list(a):
    return isinstance(a, list)


def is_dict(a):
    return isinstance(a, dict)


def is_atom(a):
    return not (hasattr(a, '__iter__'))


def get_entry(spec=None, distrib=None, distrib_version=None, pkg=None,
              section=None):
    """
    Get one entry with precedence distrib with version > distrib >
    match distrib > wildcard.
    """

    distrib_full = '{0}-{1}'.format(distrib, distrib_version)

    if pkg in spec[section]:

        if distrib_full in spec[section][pkg]:
            return spec[section][pkg][distrib_full]

        elif distrib in spec[section][pkg]:
            return spec[section][pkg][distrib]

        elif distrib in spec['match']:
            match_distrib = spec['match'][distrib]
            if match_distrib in spec[section][pkg]:
                return spec[section][pkg][match_distrib]

        if wildcard(spec) in spec[section][pkg]:
            return spec[section][pkg][wildcard(spec)]

        else:
            return None

def pkg_entries(spec=None, distrib=None, distrib_version=None, pkg=None):
    """
    Find recursively entries for pkg and distribution distrib in a
    specification spec.
    """

    result = None

    if pkg in spec['pkgs']:
        result = get_entry(spec, distrib, distrib_version, pkg, 'pkgs')

    if result is None:
        return [result]

    elif is_dict(result):
        return [result]
    else:

        if is_atom(result):
            result = [result]
        r = list()

        for e in result:

            if e != pkg:
                ne = pkg_entries(spec=spec, distrib=distrib,
                                 distrib_version=distrib_version, pkg=e)

                if ne == [None]:
                    r.append(e)

                else:
                    r += ne

        return r


def begin(distrib=None, distrib_version=None, output_mode=None):
    """
    Distribution preamble.
    """
    if output_mode == OutputMode.Docker:
        sys.stdout.write('FROM {0}:{1}\n'.format(distrib, distrib_version))


def env(definitions=None, output_mode=None):
    """
    Environment specification.
    """

    if len(definitions) > 0:

        items = list()

        if output_mode == OutputMode.Docker:
            items.append('ENV')

        items += definitions

        sys.stdout.write('{0}\n'.format(' \\ \n  '.join(items)))


def install(installer=None, command=None, pkg=None, pkgs=None,
            output_mode=OutputMode.Script):
    """
    Format an install command according to output mode.
    """

    if output_mode == OutputMode.Docker:
        items = ['RUN']

    else:
        if output_mode == OutputMode.Script:
            items = []

        else:
            sys.stderr.write('output mode {0} is not implemented\n'.format(
                output_mode))
            exit(1)

    if installer is not None:
        items.append(installer)

    if command is not None:
        if '&&' in command:
            coms = command.split('&&')
            items += ['{0} &&'.format(c.lstrip().rstrip()) for c in coms[:-1]]
            items.append(coms[-1].lstrip().rstrip())
        else:
            items.append(command)

    if pkg is not None:
        items.append(pkg)

    if pkgs is not None:
        items += pkgs

    sys.stdout.write('{0}\n'.format(' \\ \n  '.join(items)))


class Options(object):

    def __init__(self, *args, **kwargs):
        if len(args) > 0:
            for arg in args:
                self.set_options(arg)

        for k in kwargs:
            setattr(self, k, kwargs[k])

    def set_options(self, values):
        for k in values:
            setattr(self, k, values[k])

    def all(self):
        options = (name for name in dir(self) if not name.startswith('_'))
        values = (self.o for o in options)
        return dict(zip(options, values))


def get_options():

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                       ['pkg=',
                                        'pkgs=',
                                        'script',
                                        'docker',
                                        'vagrant',
                                        'split=',
                                        'distrib='])

    except getopt.GetoptError as err:
            sys.stderr.write('{0}\n'.format(str(err)))
            usage()
            exit(2)

    options = Options()

    options.distrib = None
    options.distrib_version = None
    options.pkgs = list()

    options.output_mode = OutputMode.Script
    options.split = False

    for o, a in opts:
        if o == '--distrib':
            if ':' in a:
                options.distrib, options.distrib_version = a.split(':')
            else:
                options.distrib = a
        elif o == '--pkg':
            options.pkgs.append(a)
        elif o == '--pkgs':
            options.pkgs += a.split(',')
        elif o == '--script':
            options.output_mode = OutputMode.Script
        elif o == '--docker':
            options.output_mode = OutputMode.Docker
        elif o == '--vagrant':
            options.output_mode = OutputMode.Vagrant
        elif o == '--split':
            options.split = a.lower() in ['true', 'yes', '1']

    options.specfilename = args[0]

    return options


def print_commands(*args, **kwargs):

    if len(args) == 1:
        options = args[0]
    else:
        options = Options(kwargs)

    with open(options.specfilename) as specfile:

        spec = yaml.load(specfile.read())

        by_installer = list()
        by_command = list()
        definitions = list()

        for pkg in options.pkgs:

            definition = get_entry(spec, options.distrib,
                                   options.distrib_version, pkg, 'env')

            if definition is not None:
                if is_list(definition):
                    for iter_def in definition:
                        definitions.append(iter_def)
                else:
                    definitions.append(definition)

            entries = pkg_entries(spec=spec, distrib=options.distrib,
                                  distrib_version=options.distrib_version,
                                  pkg=pkg)

            for entry in entries:
                if entry is not None:
                    if hasattr(entry, 'has_key'):
                        if 'command' in entry:
                            by_command.append(entry['command'])
                    elif hasattr(entry, 'sort'):
                        by_installer += entry
                    else:
                        by_installer.append(entry)
                else:
                    by_installer.append(pkg)

        begin(distrib=options.distrib, distrib_version=options.distrib_version,
              output_mode=options.output_mode)

        installer = get_entry(spec, options.distrib, options.distrib_version,
                              wildcard(spec), 'installer')

        assert installer is not None

        updater = get_entry(spec, options.distrib, options.distrib_version,
                            wildcard(spec), 'updater')

        if updater:
            installer = '{0} && {1}'.format(updater, installer)

        if options.split:
            for pkg in by_installer:
                format(installer=installer,
                       pkg=pkg, output_mode=options.output_mode)
        else:
            install(installer=installer,
                    pkgs=by_installer, output_mode=options.output_mode)

        for command in by_command:
            install(command=command, output_mode=options.output_mode)

        env(definitions, options.output_mode)

if __name__ == '__main__':
    print_commands(get_options())
