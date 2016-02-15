#!/usr/bin/env python

#
# Make software environment from a yml database.
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


        if wildcard in spec[section][pkg]:
            return spec[section][pkg][wildcard]

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


def run(installer=None, command=None, pkg=None, pkgs=None,
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

distrib = None
distrib_version = None
pkgs = list()

output_mode = OutputMode.Script
split = False

for o, a in opts:
    if o == '--distrib':
        if ':' in a:
            distrib, distrib_version = a.split(':')
        else:
            distrib = a
    elif o == '--pkg':
        pkgs.append(a)
    elif o == '--pkgs':
        pkgs += a.split(',')
    elif o == '--script':
        output_mode = OutputMode.Script
    elif o == '--docker':
        output_mode = OutputMode.Docker
    elif o == '--vagrant':
        output_mode = OutputMode.Vagrant
    elif o == '--split':
        split = a.lower() in ['true', 'yes', '1']


specfilename = args[0]

with open(specfilename) as specfile:

    spec = yaml.load(specfile.read())

    wildcard = None
    if 'wildcard' in spec:
        wildcard = spec['wildcard']
    else:
        wildcard = 'any'

    by_installer = list()
    by_command = list()
    definitions = list()

    for pkg in pkgs:

        definition = get_entry(spec, distrib, distrib_version, pkg, 'env')

        if definition is not None:
            if is_list(definition):
                for iter_def in definition:
                    definitions.append(iter_def)
            else:
                definitions.append(definition)

        entries = pkg_entries(spec=spec, distrib=distrib,
                              distrib_version=distrib_version, pkg=pkg)

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

    begin(distrib=distrib, distrib_version=distrib_version,
          output_mode=output_mode)

    installer = get_entry(spec, distrib, distrib_version, wildcard,
                          'installer')

    assert installer is not None

    updater = get_entry(spec, distrib, distrib_version, wildcard, 'updater')

    if updater:
        installer = '{0} && {1}'.format(updater, installer)

    if split:
        for pkg in by_installer:
            run(installer=installer,
                pkg=pkg, output_mode=output_mode)
    else:
        run(installer=installer,
            pkgs=by_installer, output_mode=output_mode)

    for command in by_command:
        run(command=command, output_mode=output_mode)

    env(definitions, output_mode)
