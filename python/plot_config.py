"""
Import this file to configure plot backend and deal with tests or remote
runtime context.

 Siconos is a program dedicated to modeling, simulation and control
 of non smooth dynamical systems.

 Copyright 2026 INRIA.

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
