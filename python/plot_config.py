# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

"""
Import this file to configure plot backend and deal with tests or remote
runtime context.

 Siconos is a program dedicated to modeling, simulation and control
 of non smooth dynamical systems.

 SPDX-FileCopyrightText: INRIA
 SPDX-License-Identifier: Apache-2.0
"""

import os
import matplotlib

"""
Test the environment and determine if plot
(e.g. plt.show()) is activated or not.

This is useful for examples that may be
called during cmake tests or remotely
"""


def choose_backend(enable_plot=True):

    if "PYTEST_CURRENT_TEST" in os.environ:
        enable_plot = False

    if os.environ.get("DISPLAY") is None:
        enable_plot = False

    if not enable_plot:
        matplotlib.use("Agg", force=True)  # non interactive backend

    print("Matplotlib backend =", matplotlib.get_backend())
    print("enable_plot =", enable_plot)
    import matplotlib.pyplot as plt

    return plt, enable_plot
