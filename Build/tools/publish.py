#!/usr/bin/env python

# documentation publication on http://siconos.gforge.inria.fr

# ./publish [-r <sha>] [-u <gforge_user>] [-s <src dir ] [-b <build dir>] \
#           [-w workdir] [-m]

# -r followed by some git sha
# -u followed by gforge login
# -s followed by src directory
# -b followed by build directory
# -d : publish devel version

# example:

# to update site with current documentation
# ./publish -u bremond -b \
#  /scratch/maurice/build/maurice/siconos/master/Release -s ~/src/git/siconos

# to update site with 3.5.x (=rev 3194) documentation
# ./publish -r3194 -u bremond [...]

# Note: some rsync error may occurs due to some files modes on remote site


import sys
import os
import shutil
import tempfile
import re
from subprocess import check_call
from getpass import getuser
from getopt import gnu_getopt, GetoptError


#
# exit only if not imported
#
def stop(n):
    import __main__ as main

    if hasattr(main, '__file__'):
        sys.exit(n)
    else:
        raise Exception('stop', n)


# a tempdir class to be used like 'with TempDir() as tmpdir:'
# i.e. the temp directory is erased ad the end of the block
class WorkDir():
    def __init__(self, prefix, tmp=False):
        self.name = None
        self.prefix = prefix
        self.tmp = tmp

    def __enter__(self):
        # just create prefix
        if not os.path.exists(self.prefix):
            os.makedirs(self.prefix)

        # The user of mkdtemp() is responsible for deleting the
        # temporary directory and its contents when done with it.
        if self.tmp:
            self.name = tempfile.mkdtemp(prefix=self.prefix)
            return self.name
        else:
            return self.prefix

    def __exit__(self, xtype, value, traceback):
        # So we remove directory here
        if self.tmp:
            shutil.rmtree(self.name)
        else:
            pass

devel = False
user = getuser()
srcdir = None
builddir = None
revision = 'HEAD'
main_doc = False

workdir_path = '/tmp/{0}/publish'.format(getuser())

try:
    opts, args = gnu_getopt(sys.argv[1:], 'r:u:s:b:w:dm', ['devel',
                                                           'main',
                                                           'revision=',
                                                           'user=',
                                                           'srcdir=',
                                                           'builddir=',
                                                           'workdir='])

    for o, a in opts:
        if o in ['-m', '--main']:
            main_doc = True

        if o in ['-r', '--revision']:
            revision = a

        elif o in ['-u', '--user']:
            user = a

        elif o in ['-d', '--devel']:
            devel = True

        elif o in ['-s', '--srcdir']:
            srcdir = a

        elif o in ['-b', '--builddir']:
            builddir = a

        elif o in ['-w', '--workdir']:
            workdir_path = a

except GetoptError, err:
    # print help information and exit:
    sys.stderr.write(str(err))  # will print something like 'option
    # -a not recognized'
    stop(2)


def get_version(path):

    with open(os.path.join(path, 'cmake', 'SiconosVersion.cmake')) as\
            cmakefile:
        cmakefile_as_str = cmakefile.read()
        majorm = re.findall(r'MAJOR_VERSION (\w+).*', cmakefile_as_str)
        minorm = re.findall(r'MINOR_VERSION (\w+).*', cmakefile_as_str)
        patchm = re.findall(r'PATCH_VERSION (\w+).*', cmakefile_as_str)
        if len(majorm) > 0:
            return '{0}.{1}.{2}'.format(majorm[0],
                                        minorm[0],
                                        patchm[0])
        else:
            return None

with WorkDir(workdir_path) as workdir:

    if builddir is None:
        builddir = os.path.join(workdir, 'build')
        try:
            os.mkdir(builddir)
        except OSError:
            pass

    if srcdir is None:
        bsrcdir = os.path.join(workdir, 'src')
        srcdir = os.path.join(workdir, 'src', 'siconos')

        try:
            os.mkdir(bsrcdir)
        except OSError:
            pass

        try:
            os.mkdir(srcdir)
        except OSError:
            pass

        # get sources
        try:
            check_call(['git', 'clone', 'git@github.com:siconos/siconos.git'],
                       cwd=bsrcdir)
        except:
            pass

    else:
        bsrcdir = os.path.dirname(srcdir)

    check_call(['git', 'checkout', revision], cwd=srcdir)

    if not devel:
        version = get_version(srcdir)
    else:
        version = 'devel'

    assert(version is not None)

    # make documentation
    check_call(['cmake', srcdir, '-DWITH_DOCUMENTATION=TRUE'] +
               ['-DWITH_{0}_DOCUMENTATION=TRUE'.format(m) for m in
                ['numerics', 'kernel', 'control', 'mechanics', 'io']],
               cwd=builddir)
    check_call(['make', 'doc'], cwd=builddir)

    # second pass for make doc
    check_call(['make', 'doc'], cwd=builddir)

    generated_doc_path = os.path.join(builddir, 'Docs', 'build', 'html')

    # change local modes
    for root, dirs, files in os.walk(generated_doc_path):
        for d in dirs:
            os.chmod(os.path.join(root, d), 0o775)
        for f in files:
            os.chmod(os.path.join(root, f), 0o664)

    os.chmod(generated_doc_path, 0o775)

    doc_path = '{0}@scm.gforge.inria.fr:/home/groups/siconos/htdocs'.\
               format(user)

    destination = os.path.join(doc_path, version)

    # upload
    check_call(['rsync', '-rlvp', '-e', 'ssh -o "StrictHostKeyChecking no"',
                generated_doc_path, destination])

    # htaccess if this is the main documentation
    if main_doc:
        htaccess_filename = os.path.join(workdir, '.htaccess')

        with open(htaccess_filename, 'w') as htaccess:
            htaccess.write(
                'redirect 301 /index.html http://siconos.gforge.inria.fr/{0}/html/index.html\n'.\
                format(version))

        check_call(['rsync', htaccess_filename, doc_path])
