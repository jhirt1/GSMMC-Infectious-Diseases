import numpy as np
import pandas as pd

import sys
import os

sys.path.append(os.path.abspath("../Generator"))
sys.path.append(os.path.abspath("../Reactions"))

import generator
import reactions


# Counters helper
def update_disease_counters(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df.loc[df["dynamic.sirvStatus"] == "I", "dynamic.infectedDays"] += 1
    df.loc[df["dynamic.sirvStatus"] == "R", "dynamic.recoveredDays"] += 1
    df.loc[df["dynamic.sirvStatus"] == "V", "dynamic.vaccinatedDays"] += 1

    return df


# SIMULATION
def run_sirv_simulation(n_runs: int, m_steps: int):

    config = generator.read_config("../Config/config_sample.yml")
    df_init = generator.intialization()

    all_runs = []
    
    for run_id in range(n_runs):

        print(f"\n=== RUN {run_id + 1}/{n_runs} ===")

        rng = np.random.default_rng(config["seed"] + run_id)

        df = df_init

        run_history = [{
                "run": run_id,
                "t": 0,
                "S": (df["dynamic.sirvStatus"] == "S").sum(),
                "I": (df["dynamic.sirvStatus"] == "I").sum(),
                "R": (df["dynamic.sirvStatus"] == "R").sum(),
                "V": (df["dynamic.sirvStatus"] == "V").sum(),
            }]

        for t in range(m_steps):
            # movement
            df = generator.jiggle_positions(
                df,
                config["rho"],
                config["sig2"],
                rng
            )

            # disease dynamic
            df = reactions.susceptible_to_infected(df, config["sig2"])
            df = reactions.infected_to_recovered(df, config["rcd"])
            df = reactions.recovered_to_susceptible(df, config["sd"])
            df = reactions.vaccinated_to_susceptible(df, config["ved"])
            df = reactions.susceptible_to_vaccinated(df, base_prob=0.005, rng=rng)

            # update counters
            df = update_disease_counters(df)

            # snapshot
            summary = {
                "run": run_id,
                "t": t+1,
                "S": (df["dynamic.sirvStatus"] == "S").sum(),
                "I": (df["dynamic.sirvStatus"] == "I").sum(),
                "R": (df["dynamic.sirvStatus"] == "R").sum(),
                "V": (df["dynamic.sirvStatus"] == "V").sum(),
            }

            run_history.append(summary)

        all_runs.extend(run_history)

    return pd.DataFrame(all_runs)