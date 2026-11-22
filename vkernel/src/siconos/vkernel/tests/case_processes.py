# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""The cases of a python test, each in its own process (a second vnative
simulation in a process crashes), run side by side: at most WORKERS at a
time, one OpenMP thread each (small scenes; many threads per process would
load the machine for nothing)."""

import concurrent.futures
import os
import subprocess

WORKERS = min(16, os.cpu_count() or 1)


def run(commands):
    """The CompletedProcess (stdout, stderr as text) of each command, in the
    order of commands."""
    env = dict(os.environ, OMP_NUM_THREADS="1")

    def one(command):
        return subprocess.run(command, capture_output=True, text=True, env=env)

    with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
        return list(pool.map(one, commands))
