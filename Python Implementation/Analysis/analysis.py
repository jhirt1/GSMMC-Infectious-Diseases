##########################################
#           Analysis / Experiments       #
##########################################
#
# Flexible parameter-sweep harness for the SIRV agent-based model.
#
# Goal (see next_steps.md): run experiments that vary how the *population* is
# initialized -- WITHOUT editing Config/config_sample.yml -- simulate each
# variant, and report how the free-growth R0 and the time-varying effective
# Re respond. Attributes we sweep:
#   - age distribution            (static.age.weights)
#   - vaccine acceptance          (static.vaccineAcceptance.probabilityTrue)
#   - comorbidity prevalence      (static.comorbidity.probabilityTrue)
#   - movement variance           (Brownian step std `sig`, hardcoded =4 in the
#                                   generator; NOT a config field -- exposed here)
#   - population density          (rho)
#
# Design: the existing Generator/Simulation entry points re-read the YAML from
# disk and run on module globals, so there is no clean injection point. Instead
# this module is a thin, self-contained DRIVER that *reuses the pipeline's pure
# functions* (the generator samplers, all five reactions, and the stats
# estimators) but feeds them an IN-MEMORY config dict. Nothing is written to
# disk, and the Brownian step `sig` becomes a first-class, tunable parameter.
#
# RUN-FROM-DIRECTORY REQUIREMENT: like the rest of the pipeline, this must be
# executed from the Analysis/ folder so the "../Config", "../Generator",
# "../Reactions" relative paths resolve.

##########################################
#         IMPORT STATEMENTS              #
##########################################

import copy
import sys
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

# Same sys.path idiom as Simulation/simulation.py so the sibling packages import.
sys.path.append(os.path.abspath("../Generator"))
sys.path.append(os.path.abspath("../Reactions"))
sys.path.append(os.path.abspath("../Simulation"))

import generator    # noqa: E402  (path hack must precede import)
import reactions    # noqa: E402
import simulation   # noqa: E402  (reused for update_disease_counters)

# NOTE: `stats`, `plotting`, and `matplotlib` are imported LAZILY inside the
# functions that need them (run_single / plot_sweep), NOT at module top level.
# Reason: parallel worker processes re-import this module on spawn, but only run
# `_simulate_one_run` (which needs generator/reactions/simulation). Keeping the
# heavy scipy/seaborn/matplotlib stack out of the import path makes worker
# startup much cheaper and improves the parallel speedup.


##########################################
#         CONFIG INJECTION               #
##########################################

# Base config path -- resolves when run from Analysis/ (sibling of Config/).
BASE_CONFIG_PATH = "../Config/config_sample.yml"

# Default Brownian movement step std. This is the value hardcoded inside
# generator.update_positions (`sig = 4`); we surface it here so it can be swept.
DEFAULT_SIG = 4.0


def _set_nested(d: dict, dotted_key: str, value) -> None:
    """
    Set a value in a nested dict using a dotted path, e.g.
    "syntheticPopulation.static.age.weights". Mutates `d` in place.
    """
    keys = dotted_key.split(".")
    node = d
    for k in keys[:-1]:
        node = node[k]
    node[keys[-1]] = value


def build_config(overrides: dict | None = None) -> dict:
    """
    Load the base YAML config and return a DEEP COPY with `overrides` applied.

    `overrides` maps dotted config paths to replacement values, e.g.
        {"rho": 1.0,
         "syntheticPopulation.static.age.weights": [0.1, 0.3, 0.3, 0.3]}

    The base file on disk is never modified -- this is the whole point of the
    harness (experiment without touching config_sample.yml).
    """
    base = generator.read_config(BASE_CONFIG_PATH)
    cfg = copy.deepcopy(base)
    for dotted_key, value in (overrides or {}).items():
        _set_nested(cfg, dotted_key, value)
    return cfg


##########################################
#         POPULATION GENERATION          #
##########################################

def generate_population(cfg: dict, seed: int | None = None) -> pd.DataFrame:
    """
    Build an initialized population DataFrame from an in-memory config.

    Reuses the generator's per-attribute samplers (which read the module-level
    `generator.rng` and `generator.infection_risk_multipliers` globals), so we
    set those globals from `cfg` first. Mirrors
    generator.generate_synthetic_population() but takes the config as an
    argument and does NOT write a CSV.
    """
    if seed is None:
        seed = cfg["seed"]

    # Wire up the generator globals its samplers depend on.
    generator.config = cfg
    generator.infection_risk_multipliers = cfg["infectionRiskMultipliers"]
    generator.set_rng_seed(seed)

    pop_size = cfg["populationSize"]
    static_cfg = cfg["syntheticPopulation"]["static"]
    age_config = static_cfg["age"]
    comorbidity_config = static_cfg["comorbidity"]
    social_activity_config = static_cfg["socialActivity"]
    geography_config = static_cfg["geography"]
    mobility_config = static_cfg["mobility"]
    vaccine_acceptance_config = static_cfg["vaccineAcceptance"]

    rows = []
    for _ in range(pop_size):
        age = generator.generate_age_str(age_config)
        comorbidity = generator.generate_comorbidity(comorbidity_config)
        social_activity = generator.generate_social_activity(social_activity_config)
        geography = generator.generate_geography(geography_config)
        mobility = generator.generate_mobility(mobility_config)
        vaccine_acceptance = generator.generate_vaccine_acceptance(vaccine_acceptance_config)

        rows.append({
            "static.guid": generator.generate_guid(),

            "static.age": age,
            "static.ageRiskMultiplier": generator.get_risk_multiplier("age", age),

            "static.comorbidity": comorbidity,
            "static.comorbidityRiskMultiplier": generator.get_risk_multiplier("comorbidity", comorbidity),

            "static.socialActivity": social_activity,
            "static.socialActivityRiskMultiplier": generator.get_risk_multiplier("socialActivity", social_activity),

            "static.geography": geography,
            "static.geographyRiskMultiplier": generator.get_risk_multiplier("geography", geography),

            "static.mobility": mobility,
            "static.mobilityRiskMultiplier": generator.get_risk_multiplier("mobility", mobility),

            "static.vaccineAcceptance": vaccine_acceptance,
            "static.vaccineAcceptanceRiskMultiplier": generator.get_risk_multiplier("vaccineAcceptance", vaccine_acceptance),

            "dynamic.vaccineStatus": False,
            "dynamic.proactiveVaccine": False,
            "dynamic.numberOfInfections": 0,
            "dynamic.sirvStatus": "S",
            "dynamic.infectedDays": 0,
            "dynamic.vaccinatedDays": 0,
            "dynamic.recoveredDays": 0,
            "dynamic.currentLocation.xcor": 0.0,
            "dynamic.currentLocation.ycor": 0.0,
        })

    df = pd.DataFrame(rows)

    # Seed the initial infected (same logic as generator.infect_initial_population).
    n_infect = cfg["initialInfected"]
    infection_idx = generator.rng.choice(df.index, size=n_infect, replace=False).tolist()
    df.loc[infection_idx, "dynamic.sirvStatus"] = "I"

    return df


##########################################
#         MOVEMENT (sig-aware)           #
##########################################

def local_jiggle(df: pd.DataFrame, rho: float, sig: float, rng) -> pd.DataFrame:
    """
    One Brownian-motion displacement with periodic boundaries.

    Port of generator.jiggle_positions / update_positions, but takes an explicit
    `sig` (the step std) instead of the hardcoded `sig = 4`, and drops the debug
    print. `sig` is the movement-variance knob we sweep. The domain side length
    L = sqrt(N / rho), so `rho` controls density / spatial mixing.
    """
    df = df.copy()

    x = df["dynamic.currentLocation.xcor"].to_numpy(dtype=float)
    y = df["dynamic.currentLocation.ycor"].to_numpy(dtype=float)

    N = len(df)
    L = np.sqrt(N / rho)

    dx = rng.normal(loc=0.0, scale=sig, size=N)
    dy = rng.normal(loc=0.0, scale=sig, size=N)

    df["dynamic.currentLocation.xcor"] = np.mod(x + dx, L)
    df["dynamic.currentLocation.ycor"] = np.mod(y + dy, L)

    return df


##########################################
#         SIMULATION DRIVER              #
##########################################

def _snapshot(df: pd.DataFrame, run_id: int, t: int) -> dict:
    status = df["dynamic.sirvStatus"]
    return {
        "run": run_id,
        "t": t,
        "S": int((status == "S").sum()),
        "I": int((status == "I").sum()),
        "R": int((status == "R").sum()),
        "V": int((status == "V").sum()),
    }


# Read-only inputs shared with worker processes. Set once per worker via the
# pool initializer so `df_init` is pickled once per worker, not once per run.
_WORKER_STATE = {}


def _init_worker(df_init: pd.DataFrame, cfg: dict, sig: float) -> None:
    """ProcessPoolExecutor initializer: stash the shared, read-only run inputs."""
    _WORKER_STATE["df_init"] = df_init
    _WORKER_STATE["cfg"] = cfg
    _WORKER_STATE["sig"] = sig


def _simulate_one_run(run_id: int) -> list[dict]:
    """
    Simulate a single independent epoch and return its per-step snapshot history.

    Each run is fully determined by its own `rng = default_rng(seed + run_id)`:
    `local_jiggle` and both stochastic reactions take this generator, so runs are
    independent across worker processes AND reproducible regardless of execution
    order or worker count. This is a top-level function so it is picklable by the
    `spawn` start method (importable as `analysis._simulate_one_run`).
    """
    df_init = _WORKER_STATE["df_init"]
    cfg = _WORKER_STATE["cfg"]
    sig = _WORKER_STATE["sig"]

    rng = np.random.default_rng(cfg["seed"] + run_id)
    df = df_init.copy()

    history = [_snapshot(df, run_id, 0)]

    for t in range(cfg["tSpan"]):
        df = local_jiggle(df, cfg["rho"], sig, rng)

        df = reactions.susceptible_to_infected(df, cfg["sig2"], rng=rng)
        df = reactions.infected_to_recovered(df, cfg["rcd"], rng=rng)
        df = reactions.recovered_to_susceptible(df, cfg["sd"], rng=rng)
        df = reactions.vaccinated_to_susceptible(df, cfg["ved"], rng=rng)
        df = reactions.susceptible_to_vaccinated(df, rng=rng)

        df = simulation.update_disease_counters(df)
        history.append(_snapshot(df, run_id, t + 1))

    return history


def run_simulation(df_init: pd.DataFrame, cfg: dict, sig: float = DEFAULT_SIG,
                   n_jobs: int | None = None) -> pd.DataFrame:
    """
    Run nSim epochs of tSpan daily steps on a fixed initial population, in parallel.

    Parameterised local equivalent of simulation.run_sirv_simulation()'s loop. The
    nSim runs are independent (each driven by its own seeded `rng`), so they are
    distributed across processes with a ProcessPoolExecutor. `n_jobs` controls the
    worker count: None -> all CPU cores; 1 -> serial (no pool overhead). Because
    every run is determined solely by `seed + run_id`, the output is identical for
    any `n_jobs` and reproducible across calls.

    Returns the RAW per-run long (run, t, S, I, R, V) DataFrame -- nSim rows per
    timestep, preserving between-run variance. Pass it through `average_runs` for
    the population-level averaged epidemic curve.
    """
    n_runs = cfg["nSim"]

    workers = (os.cpu_count() or 1) if n_jobs is None else n_jobs
    workers = max(1, min(workers, n_runs))

    if workers == 1:
        # Serial path: no process overhead, easy to debug.
        _init_worker(df_init, cfg, sig)
        histories = [_simulate_one_run(run_id) for run_id in range(n_runs)]
    else:
        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_init_worker,
            initargs=(df_init, cfg, sig),
        ) as ex:
            # ex.map preserves input order, so the result is deterministic.
            histories = list(ex.map(_simulate_one_run, range(n_runs)))

    all_runs = []
    for history in histories:
        all_runs.extend(history)

    return pd.DataFrame(all_runs)


def average_runs(sim_df: pd.DataFrame) -> pd.DataFrame:
    """
    Population-level average over runs: mean S/I/R/V per timestep across the
    nSim runs. Returns one row per `t` with columns (t, S, I, R, V).

    This is the count-level equivalent of simulation_avg.ipynb's "divide SIRV
    counts by n_runs" -- the mean over the nSim run rows at each t. Column names
    are kept as S/I/R/V so stats.py and plotting.py consume the result unchanged.
    """
    return (sim_df.groupby("t")[["S", "I", "R", "V"]]
            .mean()
            .reset_index())


##########################################
#         EXPERIMENT HARNESS             #
##########################################

# The Brownian movement std is not a config field, so we sweep it through the
# `sig` argument rather than via build_config. This sentinel marks that case.
MOVEMENT_PARAM = "movement.sig"


def run_single(overrides: dict | None = None, sig: float = DEFAULT_SIG,
               seed: int | None = None, n_jobs: int | None = None) -> dict:
    """
    Generate one population from `overrides`, simulate it, and compute metrics.

    `n_jobs` is forwarded to run_simulation (None -> all cores, 1 -> serial).

    Returns {"df", "df_avg", "cfg", "sig", "r0", "re_series", "re_peak"} where
    `df` is the raw per-run frame and `df_avg` is the population-level run
    average (mean S/I/R/V per t). Metrics are computed on `df_avg`. NOTE:
    `re_peak` is the raw max of the renewal-equation Re, which is noisy -- it is
    dominated by the sharp spike at epidemic takeoff (the first value after the
    rcd-day burn-in, divided by a near-zero prior-window mean). Treat it as
    indicative.
    """
    import stats  # lazy: keep scipy out of the parallel-worker import path

    cfg = build_config(overrides)
    df_init = generate_population(cfg, seed=seed)
    sim_df = run_simulation(df_init, cfg, sig=sig, n_jobs=n_jobs)
    df_avg = average_runs(sim_df)

    r0 = stats.calculate_r0(df_avg, rcd=cfg["rcd"])
    re_series = stats.calculate_re(df_avg, rcd=cfg["rcd"], ved=cfg["ved"])
    re_peak = float(np.nanmax(re_series["Re"].to_numpy())) if len(re_series) else np.nan

    return {
        "df": sim_df,
        "df_avg": df_avg,
        "cfg": cfg,
        "sig": sig,
        "r0": float(r0),
        "re_series": re_series,
        "re_peak": re_peak,
    }


def run_experiment(param: str, values: list, base_overrides: dict | None = None,
                   sig: float = DEFAULT_SIG, seed: int | None = None,
                   n_jobs: int | None = None) -> dict:
    """
    Sweep one attribute across `values`, returning {value: run_single(...) result}.

    `param` is either a dotted config path (applied via build_config) or the
    special MOVEMENT_PARAM sentinel ("movement.sig"), in which case each value is
    passed as the Brownian step `sig`. `base_overrides` lets you hold other
    attributes at a non-default setting while sweeping `param`. `n_jobs` is
    forwarded to each run_single's simulation (parallelism is per-run, while the
    sweep over `values` stays serial).

    NOTE: `values` for an age/social/etc. weight sweep are full weight LISTS,
    e.g. run_experiment("syntheticPopulation.static.age.weights",
                        [[0.7,0.1,0.1,0.1], [0.1,0.1,0.1,0.7]]).
    """
    results = {}
    for value in values:
        overrides = dict(base_overrides or {})
        this_sig = sig
        if param == MOVEMENT_PARAM:
            this_sig = value
        else:
            overrides[param] = value
        # Use a hashable key (lists -> tuple) so results can be dict-keyed.
        key = tuple(value) if isinstance(value, list) else value
        results[key] = run_single(overrides, sig=this_sig, seed=seed, n_jobs=n_jobs)
    return results


##########################################
#         PLOTTING                       #
##########################################

def plot_sweep(results: dict, param_label: str, plot_curves: bool = False) -> None:
    """
    Summarise a sweep: free-growth R0 and peak effective Re vs the swept value.

    If `plot_curves` is True, also overlay the SIRV trajectories for each value
    using the existing plotting.plot_sirv_onehot helper.
    """
    import matplotlib.pyplot as plt  # lazy: keep matplotlib out of worker startup
    import plotting

    keys = list(results.keys())
    # x-axis: use the value directly if scalar, else an index (e.g. weight lists).
    if all(isinstance(k, (int, float)) for k in keys):
        xs = list(keys)
        xticklabels = None
    else:
        xs = list(range(len(keys)))
        xticklabels = [str(k) for k in keys]

    r0s = [results[k]["r0"] for k in keys]
    re_peaks = [results[k]["re_peak"] for k in keys]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].plot(xs, r0s, "o-", color="tab:red")
    axes[0].set_title(f"Free-growth $R_0$ vs {param_label}")
    axes[0].set_xlabel(param_label)
    axes[0].set_ylabel("$R_0$")

    axes[1].plot(xs, re_peaks, "s-", color="tab:blue")
    axes[1].set_title(f"Peak $R_e$ vs {param_label}")
    axes[1].set_xlabel(param_label)
    axes[1].set_ylabel("peak $R_e$")

    for ax in axes:
        if xticklabels is not None:
            ax.set_xticks(xs)
            ax.set_xticklabels(xticklabels, rotation=30, ha="right")
        ax.grid(alpha=0.3)

    fig.tight_layout()
    plt.show()

    if plot_curves:
        for k in keys:
            # Use the run-averaged curve so the y-axis is on population scale.
            plotting.plot_sirv_onehot(results[k]["df_avg"], title=f"SIRV -- {param_label} = {k}")


##########################################
#         CONVENIENCE WRAPPERS           #
##########################################
# Each maps the conceptual knob to its dotted config path (or the movement
# sentinel) so notebook code reads cleanly.

def sweep_age(weight_lists: list, **kw) -> dict:
    """weight_lists: list of 4-element age weight vectors (categories 0-17,18-49,50-64,65+)."""
    return run_experiment("syntheticPopulation.static.age.weights", weight_lists, **kw)


def sweep_vaccine_acceptance(probs: list, **kw) -> dict:
    return run_experiment("syntheticPopulation.static.vaccineAcceptance.probabilityTrue", probs, **kw)


def sweep_comorbidity(probs: list, **kw) -> dict:
    return run_experiment("syntheticPopulation.static.comorbidity.probabilityTrue", probs, **kw)


def sweep_movement(sigs: list, **kw) -> dict:
    return run_experiment(MOVEMENT_PARAM, sigs, **kw)


def sweep_density(rhos: list, **kw) -> dict:
    return run_experiment("rho", rhos, **kw)


##########################################
#         DEMO ENTRY POINT               #
##########################################

if __name__ == "__main__":
    # Small, fast demo so `python analysis.py` exercises the whole path. Shrink
    # the run via overrides so this is quick; real experiments use full config.
    fast = {"nSim": 2, "tSpan": 60, "populationSize": 300}

    print("Density sweep (rho)...")
    res = sweep_density([0.25, 0.5, 1.0], base_overrides=fast)
    for k, v in res.items():
        print(f"  rho={k}: R0={v['r0']:.3f}, peak Re={v['re_peak']:.3f}")
