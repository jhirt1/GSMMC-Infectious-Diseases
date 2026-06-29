import pandas as pd
import numpy as np

def infected_to_recovered(df: pd.DataFrame, rcd: float, rng=None) -> pd.DataFrame:
    '''
    Converts an infected individual to recovered, with the infectious period drawn
    from an exponential distribution of mean `rcd` days. The individual does not
    join the susceptible pool.

    As with the other waning transitions, the exponential is memoryless, so we
    apply the equivalent per-day hazard: each infected individual recovers on a
    given day with probability p = 1 - exp(-1/rcd), yielding a geometric infectious
    period with mean `rcd` (the discrete analogue of Exp(1/rcd)).

    Pass an explicit `rng` (a numpy Generator) for reproducible / parallel-safe
    draws; if None, falls back to the global `np.random` (legacy behaviour).
    '''
    result = df.copy()

    if rng is None:
        rng = np.random

    infected_mask = (result['dynamic.sirvStatus'] == 'I').values
    n = infected_mask.sum()
    if n == 0:
        return result

    # Per-day probability of recovering for an Exp(mean=rcd) infectious period.
    p = 1.0 - np.exp(-1.0 / rcd)
    recovers = rng.random(n) < p

    newly_recovered_idx = result.index[np.where(infected_mask)[0][recovers]]
    result.loc[newly_recovered_idx, 'dynamic.sirvStatus'] = 'R'
    result.loc[newly_recovered_idx, "dynamic.recoveredDays"] = 0
    result.loc[newly_recovered_idx, "dynamic.infectedDays"] = 0

    return result

def infection_probability(df: pd.DataFrame, sig2: float) -> np.ndarray:
    '''
    Compute infection probability for each susceptible individual based on proximity and infectivity of infected individuals.
    Returns an array of probabilities aligned to the full DataFrame index.
    '''
    tsig2 = 2 * sig2

    susceptible_mask = df['dynamic.sirvStatus'] == 'S'
    infected_mask = df['dynamic.sirvStatus'] == 'I'

    susceptible_idx = df.index[susceptible_mask]
    infected_idx = df.index[infected_mask]

    pi_full = np.zeros(len(df))

    if len(susceptible_idx) == 0 or len(infected_idx) == 0:
        return pi_full

    xi = df.loc[susceptible_idx, 'dynamic.currentLocation.xcor'].values
    yi = df.loc[susceptible_idx, 'dynamic.currentLocation.ycor'].values

    xj = df.loc[infected_idx, 'dynamic.currentLocation.xcor'].values
    yj = df.loc[infected_idx, 'dynamic.currentLocation.ycor'].values
    aRM = df.loc[susceptible_idx, 'static.ageRiskMultiplier'].values
    cRM = df.loc[susceptible_idx, 'static.comorbidityRiskMultiplier'].values
    saRM = df.loc[susceptible_idx, 'static.socialActivityRiskMultiplier'].values
    gRM = df.loc[susceptible_idx, 'static.geographyRiskMultiplier'].values
    mRM = df.loc[susceptible_idx, 'static.mobilityRiskMultiplier'].values
    vaRM = df.loc[susceptible_idx, 'static.vaccineAcceptanceRiskMultiplier'].values

    # Pairwise distances: shape (n_susceptible, n_infected)
    dx = xi[:, None] - xj[None, :]
    dy = yi[:, None] - yj[None, :]
    rij2 = dx**2 + dy**2

    # Sum infectivity contributions for each susceptible individual
    pi = (aRM[:, None] * cRM[:, None] * saRM[:, None] * gRM[:, None] * mRM[:, None] * vaRM[:, None] * np.exp(-rij2 / tsig2)).sum(axis=1)

    pi_full[susceptible_mask.values] = pi
    return pi_full


def susceptible_to_infected(df: pd.DataFrame, sig2: float, rng=None) -> pd.DataFrame:
    '''
    Update infection status from susceptible to infected based on stochastic threshold applied to infection probabilities.

    Pass an explicit `rng` (a numpy Generator) for reproducible / parallel-safe
    draws; if None, falls back to the global `np.random` (legacy behaviour).
    '''
    result = df.copy()

    if rng is None:
        rng = np.random

    pi = infection_probability(result, sig2)

    susceptible_mask = result['dynamic.sirvStatus'] == 'S'
    susceptible_positions = np.where(susceptible_mask.values)[0]

    if len(susceptible_positions) == 0:
        return result

    # Stochastic threshold: infect if floor(pi + U[0,1)) >= 1
    pi_susceptible = pi[susceptible_positions]
    newly_infected = np.floor(pi_susceptible + rng.random(len(susceptible_positions))).astype(bool)

    newly_infected_idx = result.index[susceptible_positions[newly_infected]]
    result.loc[newly_infected_idx, 'dynamic.sirvStatus'] = 'I'
    result.loc[newly_infected_idx, "dynamic.infectedDays"] = 0

    return result

def recovered_to_susceptible(df: pd.DataFrame, sd: float, rng=None) -> pd.DataFrame:
    '''
    Converts a recovered individual back to susceptible, with the immunity
    duration drawn from an exponential distribution of mean `sd` days.

    The exponential is memoryless, so rather than tracking a per-person draw we
    apply the equivalent per-day hazard: each recovered individual loses immunity
    on a given day with probability p = 1 - exp(-1/sd). Over many days this yields
    a geometric waiting time with mean `sd`, the discrete analogue of Exp(1/sd).

    Pass an explicit `rng` (a numpy Generator) for reproducible / parallel-safe
    draws; if None, falls back to the global `np.random` (legacy behaviour).
    '''
    result = df.copy()

    if rng is None:
        rng = np.random

    recovered_mask = (result['dynamic.sirvStatus'] == 'R').values
    n = recovered_mask.sum()
    if n == 0:
        return result

    # Per-day probability of losing immunity for an Exp(mean=sd) dwell time.
    p = 1.0 - np.exp(-1.0 / sd)
    becomes_susceptible = rng.random(n) < p

    newly_susceptible_idx = result.index[np.where(recovered_mask)[0][becomes_susceptible]]
    result.loc[newly_susceptible_idx, 'dynamic.sirvStatus'] = 'S'
    result.loc[newly_susceptible_idx, "dynamic.recoveredDays"] = 0
    return result

# def susceptible_to_vaccinated(df: pd.DataFrame, base_prob: float = 0.005, rng=None) -> pd.DataFrame:
#     '''
#     Converts susceptible individuals to vaccinated based on demographic and behavioural multipliers.
#     '''
#     result = df.copy()

#     susceptible_mask = result['dynamic.sirvStatus'] == 'S'
#     susceptible_positions = np.where(susceptible_mask.values)[0]

#     if len(susceptible_positions) == 0:
#         return result

#     sus = result.iloc[susceptible_positions]

#     prob = (
#         base_prob
#         * sus['static.ageRiskMultiplier']
#         * sus['static.comorbidityRiskMultiplier']
#         * sus['static.socialActivityRiskMultiplier']
#         * sus['static.geographyRiskMultiplier']
#         * sus['static.mobilityRiskMultiplier']
#         * sus['static.vaccineAcceptanceRiskMultiplier']
#     ).clip(0, 1).values

#     if rng is None:
#         rng = np.random.default_rng()

#     newly_vaccinated = rng.random(len(susceptible_positions)) < prob

#     newly_vaccinated_idx = result.index[susceptible_positions[newly_vaccinated]]
#     result.loc[newly_vaccinated_idx, "dynamic.sirvStatus"] = "V"
#     result.loc[newly_vaccinated_idx, "dynamic.vaccinatedDays"] = 0

#     return result

def susceptible_to_vaccinated(df: pd.DataFrame, target_fraction: float = 0.57, n_days: int = 180, rng=None) -> pd.DataFrame:
    '''
    Converts susceptible individuals to vaccinated based on demographic and behavioural multipliers.
    '''
    result = df.copy()

    susceptible_mask = result['dynamic.sirvStatus'] == 'S'
    susceptible_positions = np.where(susceptible_mask.values)[0]

    if len(susceptible_positions) == 0:
        return result

    if rng is None:
        rng = np.random.default_rng()

    sus = result.iloc[susceptible_positions]
    w = sus['static.ageRiskMultiplier'] * sus['static.comorbidityRiskMultiplier'] * sus['static.socialActivityRiskMultiplier'] * sus['static.geographyRiskMultiplier'] * sus['static.mobilityRiskMultiplier'] * (1.0/sus['static.vaccineAcceptanceRiskMultiplier'])
    # prob = (1 - w* (1-target_fraction)**(1/n_days)).clip(0, 1).values
    prob = w*rng.poisson(0.008,len(susceptible_positions)).clip(0,1)

    newly_vaccinated = rng.random(len(susceptible_positions)) < prob

    newly_vaccinated_idx = result.index[susceptible_positions[newly_vaccinated]]
    result.loc[newly_vaccinated_idx, "dynamic.sirvStatus"] = "V"
    result.loc[newly_vaccinated_idx, "dynamic.vaccinatedDays"] = 0

    return result

def vaccinated_to_susceptible(df: pd.DataFrame, vd: float, rng=None) -> pd.DataFrame:
    '''
    Converts a vaccinated individual back to susceptible, with the vaccine-derived
    immunity duration drawn from an exponential distribution of mean `vd` days.

    As with `recovered_to_susceptible`, the exponential is memoryless, so we apply
    the equivalent per-day hazard: each vaccinated individual loses immunity on a
    given day with probability p = 1 - exp(-1/vd), yielding a geometric waiting
    time with mean `vd` (the discrete analogue of Exp(1/vd)).

    Pass an explicit `rng` (a numpy Generator) for reproducible / parallel-safe
    draws; if None, falls back to the global `np.random` (legacy behaviour).
    '''
    result = df.copy()

    if rng is None:
        rng = np.random

    vaccinated_mask = (result['dynamic.sirvStatus'] == 'V').values
    n = vaccinated_mask.sum()
    if n == 0:
        return result

    # Per-day probability of losing immunity for an Exp(mean=vd) dwell time.
    p = 1.0 - np.exp(-1.0 / vd)
    becomes_susceptible = rng.random(n) < p

    newly_susceptible_idx = result.index[np.where(vaccinated_mask)[0][becomes_susceptible]]
    result.loc[newly_susceptible_idx, 'dynamic.sirvStatus'] = 'S'
    result.loc[newly_susceptible_idx, "dynamic.vaccinatedDays"] = 0
    return result