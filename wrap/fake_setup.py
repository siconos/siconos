"""False setup.py used to find where packages installed with
'python setup.py install'
will be installed.
This is the only portable and systematic way we found.

site.getsitepackages() or distutils.sysconfig.get_python_lib() do not
behave equally on every platforms.
"""
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
# List of python modules (directories) to be included
packages = ['fake',
            ]
config = Configuration(
    name='fake',
    version='',
    packages=packages,
)

setup(**config.todict())
