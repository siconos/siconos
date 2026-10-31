# Usage: timed.py [--dataflow] [time stepping options] <script.py> [script arguments...]
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
# The time stepping options replace what the script chose (what is not given
# is left as the script set it), e.g.:
#
#   python3 timed.py --scheme at-apgd --h 1e-3 --tolerance 1e-10 bench-disks.py vnative 50
#   python3 timed.py --scheme dem --kn 1e7 --h 1e-5 bench-disks.py vnative 50
#
#   --scheme mj|gosi|at|at-apgd|dem  the scheme of the integrator
#       (simul/selectable_scheme.hpp): Moreau-Jean on the local problem, on
#       the global one, Anitescu-Tasora (projected Gauss-Seidel, APGD),
#       explicit DEM (penalty springs)
#   --theta, --gamma, --constraint-activation-threshold  of the integrator
#       (dem: theta 1 by default, the symplectic Euler scheme)
#   --h  the time step
#   --solver NAME  the solver id, a name of siconos.numerics.solver_ids,
#       SICONOS_ optional (FRICTION_2D_NSGS...); --scheme gosi needs a global
#       solver: GLOBAL_FRICTION_2D_NSN_AC if none is given
#   --tolerance, --max-iter  of the solver (also those of the
#       Anitescu-Tasora solvers)
#   --kn, --cn, --kt  the penalty law of dem: k_n 1e5, c_n 0, k_t 2/7 k_n
#       by default, unless the script gave an ExplicitDEMOSI
#
# The time stepping that ran is printed before the timings.
#
# The script is not modified: the simulations it creates time their steps,
# the bridge objects it creates take the options.

import argparse
import atexit
import runpy
import sys

import nonos.bridge as bridge
import siconos.numerics as sn

SCHEMES = {
    "mj": bridge.SCHEME_LOCAL_PROBLEM,
    "gosi": bridge.SCHEME_GLOBAL_PROBLEM,
    "at": bridge.SCHEME_CONE_COMPLEMENTARITY,
    "at-apgd": bridge.SCHEME_CONE_COMPLEMENTARITY_APGD,
    "dem": bridge.SCHEME_EXPLICIT_PENALTY,
}
SCHEME_NAMES = {v: k for k, v in SCHEMES.items()}
SOLVER_NAMES = {
    getattr(sn.solver_ids, n): n for n in dir(sn.solver_ids) if n.startswith("SICONOS_")
}

parser = argparse.ArgumentParser(
    usage="timed.py [--dataflow] [time stepping options] <script.py> [args...]",
    allow_abbrev=False,
)
parser.add_argument("--dataflow", action="store_true")
parser.add_argument("--scheme", choices=SCHEMES)
parser.add_argument("--theta", type=float)
parser.add_argument("--gamma", type=float)
parser.add_argument("--constraint-activation-threshold", type=float)
parser.add_argument("--h", type=float)
parser.add_argument("--solver")
parser.add_argument("--tolerance", type=float)
parser.add_argument("--max-iter", type=int)
parser.add_argument("--kn", type=float)
parser.add_argument("--cn", type=float)
parser.add_argument("--kt", type=float)
parser.add_argument("script")
parser.add_argument("args", nargs=argparse.REMAINDER)
options = parser.parse_args()

solver = None
if options.solver:
    name = options.solver
    name = name if name.startswith("SICONOS_") else "SICONOS_" + name
    if not hasattr(sn.solver_ids, name):
        parser.error(f"unknown solver {name}")
    solver = getattr(sn.solver_ids, name)
if options.scheme == "gosi":
    if solver is None:
        solver = sn.solver_ids.SICONOS_GLOBAL_FRICTION_2D_NSN_AC
    elif "GLOBAL" not in SOLVER_NAMES[solver]:
        parser.error(f"--scheme gosi needs a global solver, not {SOLVER_NAMES[solver]}")

# the time step: in the time discretisation, before the simulation reads it

time_discretisation_init = bridge.TimeDiscretisation.__init__


def time_discretisation(self, t0, h):
    time_discretisation_init(self, t0, h if options.h is None else options.h)


bridge.TimeDiscretisation.__init__ = time_discretisation

# the solver: its options, before the problem is created from them

osnspb_init = bridge.OSNSPB.__init__


def osnspb(self, dim, solvopts):
    if solver is not None and (solvopts is None or solvopts.solverId != solver):
        solvopts = sn.SolverOptions(solver)
    osnspb_init(self, dim, solvopts)
    if options.tolerance is not None:
        self._so.set_dparam(sn.params.SICONOS_DPARAM_TOL, options.tolerance)
    if options.max_iter is not None:
        self._so.set_iparam(sn.params.SICONOS_IPARAM_MAX_ITER, options.max_iter)


bridge.OSNSPB.__init__ = osnspb

# the integrator: the options over what the script chose


def set_integrator(simulation, osi):
    vosi = simulation.handle().one_step_integrator()
    if options.scheme is not None:
        vosi.set_scheme(SCHEMES[options.scheme])
    if options.scheme == "dem" and type(osi).__name__ != "ExplicitDEMOSI":
        # the defaults of bridge.ExplicitDEMOSI
        kn = 1e5 if options.kn is None else options.kn
        vosi.set_theta(1.0)
        vosi.set_normal_stiffness(kn)
        vosi.set_normal_damping(0.0)
        vosi.set_tangential_stiffness(2.0 / 7.0 * kn)
    for value, set_value in (
        (options.theta, vosi.set_theta),
        (options.gamma, vosi.set_gamma),
        (
            options.constraint_activation_threshold,
            vosi.set_constraint_activation_threshold,
        ),
        (options.kn, vosi.set_normal_stiffness),
        (options.cn, vosi.set_normal_damping),
        (options.kt, vosi.set_tangential_stiffness),
    ):
        if value is not None:
            set_value(value)


insert_integrator = bridge.Simulation.insertIntegrator


def timed_insert_integrator(self, osi):
    insert_integrator(self, osi)
    set_integrator(self, osi)


bridge.Simulation.insertIntegrator = timed_insert_integrator


def time_stepping(simulation):
    """the time stepping that ran"""
    vosi = simulation.handle().one_step_integrator()
    scheme = SCHEME_NAMES.get(vosi.scheme(), str(vosi.scheme()))
    lines = [
        f"scheme {scheme}, h {simulation._timedisc.handle().h()}, "
        f"theta {vosi.theta()}, gamma {vosi.gamma()}, "
        f"constraint activation threshold {vosi.constraint_activation_threshold()}"
    ]
    if scheme == "dem":
        lines.append(
            f"k_n {vosi.normal_stiffness()}, c_n {vosi.normal_damping()}, "
            f"k_t {vosi.tangential_stiffness()}"
        )
    osnspb = getattr(simulation, "_osnspb", None)
    if osnspb is not None and scheme != "dem":
        so = osnspb._so
        # Anitescu-Tasora: its own solver, the tolerance and the maximum
        # number of iterations from the options
        name = {"at": "projected Gauss-Seidel", "at-apgd": "APGD"}.get(
            scheme, SOLVER_NAMES.get(so.solver_id(), so.solver_id())
        )
        lines.append(
            f"solver {name}, "
            f"tolerance {so.dparam(sn.params.SICONOS_DPARAM_TOL)}, "
            f"max iterations {so.iparam(sn.params.SICONOS_IPARAM_MAX_ITER)}"
        )
    return "time stepping: " + "\n               ".join(lines)


simulation_init = bridge.Simulation.__init__


def timed_init(self, *args, **kwargs):
    simulation_init(self, *args, **kwargs)
    # without an integrator inserted by the script, the options all the same
    set_integrator(self, None)
    handle = self.handle()
    handle.enable_timings(True)
    handle.enable_dataflow(options.dataflow)

    def report():
        if options.dataflow:
            print("\n" + handle.dataflow_report(), file=sys.stderr)
        print("\n" + time_stepping(self), file=sys.stderr)
        print(handle.timings_report(), file=sys.stderr)

    atexit.register(report)


bridge.Simulation.__init__ = timed_init

sys.argv = [options.script] + options.args
runpy.run_path(options.script, run_name="__main__")
