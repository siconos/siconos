"""
Tools used during configuration and build process
(most of them used in CMake files, during build or runtime).

This file is to be copied into CMAKE_BINARY_DIR/share using configure_file

 Siconos is a program dedicated to modeling, simulation and control
 of non smooth dynamical systems.

 Copyright 2018 INRIA.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
"""


def parse_cmake_list(var):
    """Transform cmake list-like variables
    into python lists.

    Parameters
    ----------
    var : string
        like "a;b;c"

    Returns python list

    If var is already a list, does nothing.

    Example::

        a = parse_cmake_list("var1;var2;var3;")
        # --> a = ['var', 'var2', 'var3']

    """
    if isinstance(var, list):
        return var
    if var != "":
        res = list(set(var.split(';')))
        # list/set stuff to avoid duplicates
        # remove empty strings to avoid '-I -I' things leading to bugs
        if res.count(''):
            res.remove('')
        return res

    return []
