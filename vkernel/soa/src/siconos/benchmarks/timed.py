# Usage: timed.py [--dataflow] <script.py> [script arguments...]
#
# Runs a vnative script (bench-disks.py, bench-polytopes.py,
# tests/bench_fem_disks.py...) with the timings of the stages of its time
# steps, and prints them at the end, on stderr, one line per stage of the
# schedule of the time step (see storage/schedule.hpp and
# time_stepping::enable_timings), e.g.:
#
#   python3 timed.py bench-disks.py vnative 50
#
# With --dataflow, it also prints what the stages do: the kernels each stage
# ran, their reads and writes, a placement (device, host, no kernel) and the
# candidates for a fusion (see storage::schedule::dataflow).
#
# The script is not modified: the simulations it creates time their steps.

import atexit
import runpy
import sys

import nonos.bridge as bridge

dataflow = len(sys.argv) > 1 and sys.argv[1] == "--dataflow"
if dataflow:
    del sys.argv[1]

simulation_init = bridge.Simulation.__init__


def timed_init(self, *args, **kwargs):
    simulation_init(self, *args, **kwargs)
    handle = self.handle()
    handle.enable_timings(True)
    handle.enable_dataflow(dataflow)

    def report():
        if dataflow:
            print("\n" + handle.dataflow_report(), file=sys.stderr)
        print("\n" + handle.timings_report(), file=sys.stderr)

    atexit.register(report)


bridge.Simulation.__init__ = timed_init

script = sys.argv[1]
sys.argv = sys.argv[1:]
runpy.run_path(script, run_name="__main__")
