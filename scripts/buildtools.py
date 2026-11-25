# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

"""
Tools used during configuration and build process
(most of them used in CMake files, during build or runtime).

This file is to be copied into CMAKE_BINARY_DIR/share using configure_file

 Siconos is a program dedicated to modeling, simulation and control
 of non smooth dynamical systems.

 SPDX-FileCopyrightText: INRIA
 SPDX-License-Identifier: Apache-2.0
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
        res = list(set(var.split(";")))
        # list/set stuff to avoid duplicates
        # remove empty strings to avoid '-I -I' things leading to bugs
        if res.count(""):
            res.remove("")
        return res

    return []
