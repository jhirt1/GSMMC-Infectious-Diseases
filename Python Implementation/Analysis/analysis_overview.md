# `analysis.py` — Parameter-Sweep Experiment Harness

A walkthrough of what [analysis.py](analysis.py) does, how it builds on the
Generator and Simulation stages, and what it is meant to produce.

---

## 1. Purpose

The team pipeline (Generator → Simulation → Analysis) can simulate **one**
population — the one fixed by [`Config/config_sample.yml`](../Config/config_sample.yml).
`analysis.py` answers the question in [next_steps.md](next_steps.md): *how do
disease dynamics change when we change the population?*

It runs **controlled experiments** that vary one initialization attribute at a
time, simulate each variant, and report how two transmission metrics respond:

| Metric | Meaning | Source |
|---|---|---|
| **R₀** | Free-growth reproduction number (fully-susceptible phase) | [`stats.calculate_r0`](stats.py) (Euler–Lotka) |
| **peak Rₑ** | Peak of the time-varying effective reproduction number | [`stats.calculate_re`](stats.py) (Cori renewal) |

The five attributes it can sweep:

| Knob | Where it lives | How `analysis.py` varies it |
|---|---|---|
| Age distribution | `syntheticPopulation.static.age.weights` | in-memory config override |
| Vaccine acceptance | `...vaccineAcceptance.probabilityTrue` | in-memory config override |
| Comorbidity prevalence | `...comorbidity.probabilityTrue` | in-memory config override |
| Population density | `rho` | in-memory config override |
| **Movement variance** | Brownian step std `sig` (hardcoded `=4`, **not a config field**) | explicit function argument |

---

## 2. The core design problem (and the solution)

The existing entry points cannot be reused directly:

- [`generator.intialization()`](../Generator/generator.py) and
  [`simulation.run_sirv_simulation()`](../Simulation/simulation.py) **re-read the
  YAML from disk** every call and operate on **module-level globals**. There is
  no argument you can pass to inject a different population.
- The only way to change the population through those entry points is to **edit
  `config_sample.yml`** — which `next_steps.md` explicitly forbids.
- **Movement variance isn't even in the config.** `generator.update_positions`
  hardcodes `sig = 4` (the line `sig = np.sqrt(sig2)` is commented out). `sig2`
  is a *different* quantity — the width of the infection kernel used in
  `reactions.infection_probability` — so you cannot reach movement variance via
  config at all.

**Solution.** `analysis.py` is a thin, self-contained **driver** that *reuses the
pipeline's pure functions* (the generator samplers, all five reactions, and the
stats estimators) but feeds them an **in-memory config dictionary**. Nothing is
written to disk, the base config is never touched, and the Brownian step `sig`
becomes a first-class, tunable argument.

```
build_config(overrides)  ─►  generate_population(cfg)  ─►  run_simulation(df, cfg, sig)  ─►  stats.calculate_r0 / calculate_re
   (deep-copy YAML,            (reuse generator                (local copy of the              (R0, Re from the
    apply overrides)            samplers, no CSV)                daily reaction loop)            long SIRV frame)
```

---

## 3. How it builds on Generation

### `build_config(overrides)` — in-memory config injection
Loads the base YAML via `generator.read_config`, makes a **deep copy**, and
applies overrides addressed by **dotted paths** (e.g.
`"syntheticPopulation.static.age.weights"`). The helper `_set_nested` walks the
nested dict and sets the leaf. The on-disk file is never modified — this is what
makes "experiment without changing the config" possible.

### `generate_population(cfg, seed)` — reuse the samplers, skip the side effects
This mirrors [`generator.generate_synthetic_population()`](../Generator/generator.py)
but is parameterised and side-effect-free:

1. **Wires up the generator globals** the samplers depend on — sets
   `generator.config`, `generator.infection_risk_multipliers`, and seeds the RNG
   via `generator.set_rng_seed(seed)`. (The generator's samplers read these
   module globals rather than taking them as arguments, so we set them before
   calling.)
2. **Reuses the per-attribute samplers** unchanged —
   `generate_age_str`, `generate_comorbidity`, `generate_social_activity`,
   `generate_geography`, `generate_mobility`, `generate_vaccine_acceptance`,
   plus `get_risk_multiplier` and `generate_guid`. Each person becomes one row of
   dotted-path columns (`static.age`, `static.ageRiskMultiplier`,
   `dynamic.sirvStatus`, …), exactly matching the schema the rest of the pipeline
   expects.
3. **Seeds the initial infected** the same way as
   `generator.infect_initial_population` (random `initialInfected` agents set to
   `"I"`).
4. **Does NOT write a CSV** — unlike the global generator, which writes
   `synthetic_population2.csv`. The harness keeps everything in memory.

The result is a `df_init` DataFrame: the starting population for a run.

> **Note on positions.** New agents start at `(0, 0)`; the first `local_jiggle`
> inside `run_simulation` scatters them. (This matches the pipeline's own
> behaviour, where initial positions are set by a jiggle.)

---

## 4. How it builds on Simulation

### `local_jiggle(df, rho, sig, rng)` — the movement knob
A port of `generator.jiggle_positions` + `update_positions`, with two changes:

- It takes an **explicit `sig`** (the Brownian step standard deviation) instead
  of the hardcoded `4`. **This is the movement-variance knob.**
- It drops the per-step debug `print`.

The domain side length is `L = sqrt(N / rho)`, so **`rho` controls density** —
smaller `L` packs agents closer, raising contact rates. Positions wrap with
periodic boundaries (`np.mod`).

### `run_simulation(df_init, cfg, sig)` — the daily loop
A local copy of the loop in `simulation.run_sirv_simulation()`, but parameterised
by the in-memory `cfg` and `sig` instead of re-reading the YAML. For each of
`cfg["nSim"]` runs it reseeds an RNG (`seed + run_id`), copies `df_init`, then for
each of `cfg["tSpan"]` daily steps it:

1. `local_jiggle` — move everyone one Brownian step.
2. The **five reactions in order**, unchanged from
   [`reactions.py`](../Reactions/reactions.py):
   - `susceptible_to_infected(df, cfg["sig2"])` — distance-weighted Gaussian
     infection kernel.
   - `infected_to_recovered(df, cfg["rcd"])` — deterministic recovery after `rcd`
     days.
   - `recovered_to_susceptible(df, cfg["sd"])` — waning immunity after `sd` days.
   - `vaccinated_to_susceptible(df, cfg["ved"])` — vaccine waning after `ved`
     days.
   - `susceptible_to_vaccinated(df, rng=rng)` — stochastic vaccination.
3. `simulation.update_disease_counters(df)` — increment the I/R/V day counters.
   This helper is **imported from and reused as-is** from
   [`Simulation/simulation.py`](../Simulation/simulation.py#L15) — it is a pure,
   side-effect-free function, so there is no reason to duplicate it.
4. `_snapshot` — record `(run, t, S, I, R, V)` for that step.

It returns the **long DataFrame** with columns `run, t, S, I, R, V` — the exact
shape [`stats.py`](stats.py) and [`plotting.py`](plotting.py) already consume, so
those modules work without modification.

> **Reused, not reimplemented.** The reactions, the stats estimators, and
> `simulation.update_disease_counters` are the originals. `analysis.py` only
> re-creates the *orchestration that can't be reused* — the config plumbing, the
> population builder, and the daily loop (`run_simulation`), which is local
> because `simulation.run_sirv_simulation()` has no injection point.

---

## 5. The experiment layer

### `run_single(overrides, sig, seed)`
The atomic experiment: `build_config` → `generate_population` → `run_simulation`
→ compute metrics. Returns a dict:

```python
{"df": <long SIRV frame>, "cfg": <resolved config>, "sig": <movement std>,
 "r0": <float>, "re_series": <DataFrame t/incidence/Re>, "re_peak": <float>}
```

### `run_experiment(param, values, base_overrides, sig, seed)`
Sweeps one attribute across `values`, calling `run_single` per value and
collecting `{value: result}`. Two cases:

- **Config attribute** (e.g. `"rho"`): each value is applied as an override.
- **Movement** (`param == MOVEMENT_PARAM`, the sentinel `"movement.sig"`): each
  value is passed as `sig`, since movement isn't a config field.

`base_overrides` lets you **hold other attributes fixed at a non-default setting**
while sweeping `param` (also handy for shrinking `nSim`/`tSpan`/`populationSize`
during quick tests). For weight-vector sweeps (age), each `value` is a full list
and is stored under a tuple key.

### Convenience wrappers
Thin, readable entry points so notebook code reads well:

```python
sweep_age([...])                 # age weight vectors
sweep_vaccine_acceptance([...])  # probabilityTrue values
sweep_comorbidity([...])         # probabilityTrue values
sweep_movement([...])            # Brownian sig values
sweep_density([...])             # rho values
```

---

## 6. What it produces

### `plot_sweep(results, param_label, plot_curves=False)`
Two summary panels across the swept values:

- **Free-growth R₀ vs the swept value** (red).
- **Peak Rₑ vs the swept value** (blue).

With `plot_curves=True` it also overlays the **SIRV trajectories** for each value
via `plotting.plot_sirv_onehot`. Scalar sweeps (rho, sig, probabilities) plot on a
numeric x-axis; list sweeps (age weights) plot on indexed categorical ticks.

### Programmatic output
Every `run_*` returns plain dicts/DataFrames, so a notebook can post-process:
the long SIRV frame per run, the resolved config, R₀, and the full `Re(t)` series.

### Demo entry point
`python analysis.py` (run **from `Analysis/`**) executes a small, fast density
sweep and prints R₀ / peak Rₑ per value — a smoke test of the whole path.

---

## 7. How to run

**Run from the `Analysis/` directory** so the `../Config`, `../Generator`,
`../Reactions` relative paths resolve (same requirement as the rest of the
pipeline). In a notebook:

```python
import analysis

# Sweep density, holding the run small for speed:
results = analysis.sweep_density(
    [0.25, 0.5, 1.0],
    base_overrides={"nSim": 20, "tSpan": 150, "populationSize": 500},
)
analysis.plot_sweep(results, "population density (rho)", plot_curves=True)

# Sweep movement variance:
mv = analysis.sweep_movement([1.0, 4.0, 10.0],
                             base_overrides={"nSim": 20, "tSpan": 150})
analysis.plot_sweep(mv, "movement std (sig)")
```

Expected directions (sanity checks): denser populations (higher `rho`) and more
mixing (higher `sig`) raise **R₀**; higher vaccine acceptance lowers sustained
transmission.

---

## 8. Caveats inherited from the pipeline

`analysis.py` reuses the existing reactions and stats as-is, so it inherits their
known rough edges:

- **`peak Rₑ` is noisy.** The renewal-equation Rₑ has a sharp spike at epidemic
  takeoff (the first value after the `rcd`-day estimator burn-in, divided by a
  near-zero prior-window mean). The raw max is therefore dominated by that
  transient — treat it as **indicative**, and lean on R₀ and the SIRV curves for
  firm conclusions.
- **Partial reproducibility.** `susceptible_to_infected` and
  `susceptible_to_vaccinated` draw from the **global** `np.random` rather than the
  seeded `rng`, so those specific transitions are not reproducible from `seed`
  even though generation and movement are.
- **One authorized upstream fix.** `reactions.infection_probability` was corrected
  so the risk multipliers broadcast over the **susceptible** target
  (`[:, None]`) — required for the kernel to run at all after the multipliers were
  switched to index susceptibles. This is the only change outside `Analysis/`.
