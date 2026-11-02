"""Simulation options, applied to a driver script that runs in-process.

The `siconos` command loads this module before the driver runs when given
simulation options: it overrides what the driver script chose, before the
driver runs (the driver is run with runpy after apply() has installed its
patches). Nothing is given, nothing is changed.

Two groups of options:

  - the run options of siconos.io.mechanics_run.MechanicsHdf5Runner, shared
    by all the backends: t0, T, h, output_format. They are forced in
    run_options before run_initialize reads them.
  - the objects of nonos.bridge, the vnative backend only: the scheme of the
    integrator (MoreauJean, MoreauJeanGOSI, AnitescuTasora,
    AnitescuTasoraAPGD, DEM) and its theta, gamma and
    constraint activation threshold, the solver and its tolerance and maximum
    number of iterations, the penalty law of the explicit DEM scheme (dem_kn,
    dem_cn, dem_kt), and the observers: the timings of the stages of the time
    steps (timings) and the dataflow of the schedule (dataflow), reported on
    stderr at exit.

The vnative options import nonos.bridge lazily: they are only available when
the vnative backend is installed, and applying them on another backend
errors.

This module is loaded from its file path (not by importing the nonos
package), so that the run options above work without the compiled _nonos
module; the vnative options still import nonos.bridge the usual way.
"""

import atexit
import sys

# the schemes of the integrator (simul/selectable_scheme.hpp), mirrored from
# nonos.bridge
SCHEME_LOCAL_PROBLEM = 0
SCHEME_GLOBAL_PROBLEM = 1
SCHEME_CONE_COMPLEMENTARITY = 2
SCHEME_CONE_COMPLEMENTARITY_APGD = 3
SCHEME_EXPLICIT_PENALTY = 4

SCHEMES = {
    "MoreauJean": SCHEME_LOCAL_PROBLEM,
    "MoreauJeanGOSI": SCHEME_GLOBAL_PROBLEM,
    "AnitescuTasora": SCHEME_CONE_COMPLEMENTARITY,
    "AnitescuTasoraAPGD": SCHEME_CONE_COMPLEMENTARITY_APGD,
    "DEM": SCHEME_EXPLICIT_PENALTY,
}
# the short names of the first version, still accepted
SCHEME_ALIASES = {
    "mj": "MoreauJean",
    "gosi": "MoreauJeanGOSI",
    "at": "AnitescuTasora",
    "at-apgd": "AnitescuTasoraAPGD",
    "dem": "DEM",
}
SCHEME_NAMES = {v: k for k, v in SCHEMES.items()}

# the run options of MechanicsHdf5Runner, the backends all
CROSS_BACKEND_OPTIONS = ("t0", "T", "h", "output_format")
# the objects of nonos.bridge, vnative only
VNATIVE_OPTIONS = (
    "scheme",
    "theta",
    "gamma",
    "constraint_activation_threshold",
    "solver",
    "tolerance",
    "max_iter",
    "dem_kn",
    "dem_cn",
    "dem_kt",
    "timings",
    "dataflow",
)


def _error(message):
    print(f"siconos: {message}", file=sys.stderr)
    sys.exit(2)


_solver_names_cache = None
_solver_ids_cache = None


def _solver_names():
    """the solver ids of siconos.numerics as {id: name}, built once"""
    global _solver_names_cache
    if _solver_names_cache is None:
        import siconos.numerics as sn

        _solver_names_cache = {
            getattr(sn.solver_ids, n): n
            for n in dir(sn.solver_ids)
            if n.startswith("SICONOS_")
        }
    return _solver_names_cache


def _solver_id(name):
    """the id of a solver of siconos.numerics, SICONOS_ optional in name;
    exits with a clean error on an unknown name"""
    global _solver_ids_cache
    if _solver_ids_cache is None:
        _solver_ids_cache = {v: k for k, v in _solver_names().items()}
    if not name.startswith("SICONOS_"):
        name = "SICONOS_" + name
    if name not in _solver_ids_cache:
        _error(f"unknown solver {name}")
    return _solver_ids_cache[name]


def apply(options):
    """Override the simulation parameters of the driver that follows.

    options: a dict, the keys of CROSS_BACKEND_OPTIONS (any backend) and of
    VNATIVE_OPTIONS (the vnative backend only); absent or None keys are left
    as the driver set them.
    """

    def set(value):
        # None and False (the timings/dataflow flags not given) are not set
        return value is not None and value is not False

    cross = {key: options.get(key) for key in CROSS_BACKEND_OPTIONS}
    if any(set(value) for value in cross.values()):
        _patch_run_options(cross)

    vnative = {key: options.get(key) for key in VNATIVE_OPTIONS}
    if vnative["scheme"] is not None:
        # the short names of the first version are still accepted
        vnative["scheme"] = SCHEME_ALIASES.get(vnative["scheme"], vnative["scheme"])
        if vnative["scheme"] not in SCHEMES:
            _error(f"unknown scheme {vnative['scheme']}, one of " + ", ".join(SCHEMES))
    if any(set(value) for value in vnative.values()):
        _patch_vnative(vnative)


def _patch_run_options(cross):
    # the run options: in run_initialize, before the runner reads them
    from siconos.io.mechanics_run import MechanicsHdf5Runner

    run_initialize = MechanicsHdf5Runner.run_initialize

    def run_options_run_initialize(self):
        for key, value in cross.items():
            if value is not None:
                self._run_options[key] = value
        return run_initialize(self)

    MechanicsHdf5Runner.run_initialize = run_options_run_initialize


def _patch_vnative(v):
    try:
        import nonos.bridge as bridge
    except ImportError as e:
        _error(f"the vnative simulation options need the vnative backend (nonos): {e}")
    import siconos.numerics as sn

    solver = None
    if v["solver"] is not None:
        solver = _solver_id(v["solver"])
    if v["scheme"] == "MoreauJeanGOSI":
        if solver is None:
            solver = sn.solver_ids.SICONOS_GLOBAL_FRICTION_2D_NSN_AC
        elif "GLOBAL" not in _solver_names()[solver]:
            _error(
                f"--scheme MoreauJeanGOSI needs a global solver, "
                f"not {_solver_names()[solver]}"
            )

    # the solver: its options, before the problem is created from them
    osnspb_init = bridge.OSNSPB.__init__

    def osnspb(self, dim, solvopts):
        if solver is not None and (solvopts is None or solvopts.solverId != solver):
            solvopts = sn.SolverOptions(solver)
        osnspb_init(self, dim, solvopts)
        if v["tolerance"] is not None:
            self._so.set_dparam(sn.params.SICONOS_DPARAM_TOL, v["tolerance"])
        if v["max_iter"] is not None:
            self._so.set_iparam(sn.params.SICONOS_IPARAM_MAX_ITER, v["max_iter"])

    bridge.OSNSPB.__init__ = osnspb

    # the integrator: the options over what the script chose
    def set_integrator(simulation, osi):
        vosi = simulation.handle().one_step_integrator()
        if v["scheme"] is not None:
            vosi.set_scheme(SCHEMES[v["scheme"]])
        if v["scheme"] == "DEM" and type(osi).__name__ != "ExplicitDEMOSI":
            # the defaults of bridge.ExplicitDEMOSI
            kn = 1e5 if v["dem_kn"] is None else v["dem_kn"]
            vosi.set_theta(1.0)
            vosi.set_normal_stiffness(kn)
            vosi.set_normal_damping(0.0)
            vosi.set_tangential_stiffness(2.0 / 7.0 * kn)
        for value, set_value in (
            (v["theta"], vosi.set_theta),
            (v["gamma"], vosi.set_gamma),
            (
                v["constraint_activation_threshold"],
                vosi.set_constraint_activation_threshold,
            ),
            (v["dem_kn"], vosi.set_normal_stiffness),
            (v["dem_cn"], vosi.set_normal_damping),
            (v["dem_kt"], vosi.set_tangential_stiffness),
        ):
            if value is not None:
                set_value(value)

    insert_integrator = bridge.Simulation.insertIntegrator

    def vnative_insert_integrator(self, osi):
        insert_integrator(self, osi)
        set_integrator(self, osi)

    bridge.Simulation.insertIntegrator = vnative_insert_integrator

    # the observers: the time steps time themselves, the report at exit
    if v["timings"] or v["dataflow"]:
        simulation_init = bridge.Simulation.__init__

        def observed_init(self, *args, **kwargs):
            simulation_init(self, *args, **kwargs)
            # without an integrator inserted by the script, the options all the same
            set_integrator(self, None)
            handle = self.handle()
            if v["timings"]:
                handle.enable_timings(True)
            if v["dataflow"]:
                handle.enable_dataflow(True)

            def report():
                if v["dataflow"]:
                    print("\n" + handle.dataflow_report(), file=sys.stderr)
                print("\n" + _time_stepping(self), file=sys.stderr)
                if v["timings"]:
                    print(handle.timings_report(), file=sys.stderr)

            atexit.register(report)

        bridge.Simulation.__init__ = observed_init


def _time_stepping(simulation):
    """the time stepping that ran"""
    import siconos.numerics as sn

    vosi = simulation.handle().one_step_integrator()
    scheme = SCHEME_NAMES.get(vosi.scheme(), str(vosi.scheme()))
    lines = [
        f"scheme {scheme}, h {simulation._timedisc.handle().h()}, "
        f"theta {vosi.theta()}, gamma {vosi.gamma()}, "
        f"constraint activation threshold {vosi.constraint_activation_threshold()}"
    ]
    if scheme == "DEM":
        lines.append(
            f"k_n {vosi.normal_stiffness()}, c_n {vosi.normal_damping()}, "
            f"k_t {vosi.tangential_stiffness()}"
        )
    osnspb = getattr(simulation, "_osnspb", None)
    if osnspb is not None and scheme != "DEM":
        so = osnspb._so
        # Anitescu-Tasora: its own solver, the tolerance and the maximum
        # number of iterations from the options
        name = {
            "AnitescuTasora": "projected Gauss-Seidel",
            "AnitescuTasoraAPGD": "APGD",
        }.get(scheme, _solver_names().get(so.solver_id(), so.solver_id()))
        lines.append(
            f"solver {name}, "
            f"tolerance {so.dparam(sn.params.SICONOS_DPARAM_TOL)}, "
            f"max iterations {so.iparam(sn.params.SICONOS_IPARAM_MAX_ITER)}"
        )
    return "time stepping: " + "\n               ".join(lines)
