import pandas as pd
import numpy as np

def infected_to_recovered(df: pd.DataFrame, rcd: int) -> pd.DataFrame:
    '''
    Converts an infected individual to recovered after rcd days of infection. 
    The individual does not join the susceptible pool
    '''
    result = df.copy()
    
    # Find infected individuals who have reached recovery day
    infected_mask = result['dynamic.sirvStatus'] == 'I'
    recovered_mask = result['dynamic.infectedDays'] == rcd
    
    # Update status to 'recovered' for those who've completed rcd days
    result.loc[infected_mask & recovered_mask, 'dynamic.sirvStatus'] = 'R'
    result.loc[infected_mask & recovered_mask, "dynamic.recoveredDays"] = 0
    result.loc[infected_mask & recovered_mask, "dynamic.infectedDays"] = 0
    
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
    aRM = df.loc[infected_idx, 'static.ageRiskMultiplier'].values
    cRM = df.loc[infected_idx, 'static.comorbidityRiskMultiplier'].values
    saRM = df.loc[infected_idx, 'static.socialActivityRiskMultiplier'].values
    gRM = df.loc[infected_idx, 'static.geographyRiskMultiplier'].values
    mRM = df.loc[infected_idx, 'static.mobilityRiskMultiplier'].values
    vaRM = df.loc[infected_idx, 'static.vaccineAcceptanceRiskMultiplier'].values

    # Pairwise distances: shape (n_susceptible, n_infected)
    dx = xi[:, None] - xj[None, :]
    dy = yi[:, None] - yj[None, :]
    rij2 = dx**2 + dy**2

    # Sum infectivity contributions for each susceptible individual
    pi = (aRM[None, :] * cRM[None, :] * saRM[None, :] * gRM[None, :] * mRM[None, :] * vaRM[None, :] * np.exp(-rij2 / tsig2)).sum(axis=1)

    pi_full[susceptible_mask.values] = pi
    return pi_full


def susceptible_to_infected(df: pd.DataFrame, sig2: float) -> pd.DataFrame:
    '''
    Update infection status from susceptible to infected based on stochastic threshold applied to infection probabilities.
    '''
    result = df.copy()

    pi = infection_probability(result, sig2)

    susceptible_mask = result['dynamic.sirvStatus'] == 'S'
    susceptible_positions = np.where(susceptible_mask.values)[0]

    if len(susceptible_positions) == 0:
        return result

    # Stochastic threshold: infect if floor(pi + U[0,1)) >= 1
    pi_susceptible = pi[susceptible_positions]
    newly_infected = np.floor(pi_susceptible + np.random.rand(len(susceptible_positions))).astype(bool)

    newly_infected_idx = result.index[susceptible_positions[newly_infected]]
    result.loc[newly_infected_idx, 'dynamic.sirvStatus'] = 'I'
    result.loc[newly_infected_idx, "dynamic.infectedDays"] = 0

    return result

def recovered_to_susceptible(df: pd.DataFrame, sd: int) -> pd.DataFrame:
    '''
    Converts a recovered individual to susceptible after sd days. 
    '''
    result = df.copy()

    # Find recovered individuals who have reached susceptibility day
    recovered_mask = result['dynamic.sirvStatus'] == 'R'
    susceptible_mask = result['dynamic.recoveredDays'] == sd

    # Update status to 'susceptible' for those who've completed sd days
    result.loc[recovered_mask & susceptible_mask, 'dynamic.sirvStatus'] = 'S'
    result.loc[recovered_mask & susceptible_mask, "dynamic.recoveredDays"] = 0
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

    sus = result.iloc[susceptible_positions]

    prob = (1 -
        sus['static.ageRiskMultiplier']
        * sus['static.comorbidityRiskMultiplier']
        * sus['static.socialActivityRiskMultiplier']
        * sus['static.geographyRiskMultiplier']
        * sus['static.mobilityRiskMultiplier']
        * sus['static.vaccineAcceptanceRiskMultiplier']
        * (1-target_fraction)**(1/n_days)
    ).clip(0, 1).values

    if rng is None:
        rng = np.random.default_rng()

    newly_vaccinated = rng.random(len(susceptible_positions)) < prob

    newly_vaccinated_idx = result.index[susceptible_positions[newly_vaccinated]]
    result.loc[newly_vaccinated_idx, "dynamic.sirvStatus"] = "V"
    result.loc[newly_vaccinated_idx, "dynamic.vaccinatedDays"] = 0

    return result

def vaccinated_to_susceptible(df: pd.DataFrame, vd: int) -> pd.DataFrame:
    '''
    Converts a vaccinated individual to susceptible after vd days of vaccination. 
    '''
    result = df.copy()

    # Find vaccinated individuals who have reached susceptibility day
    vaccinated_mask = result['dynamic.sirvStatus'] == 'V'
    susceptible_mask = result['dynamic.vaccinatedDays'] == vd

    # Update status to 'susceptible' for those who've completed sd days
    result.loc[vaccinated_mask & susceptible_mask, 'dynamic.sirvStatus'] = 'S'
    result.loc[vaccinated_mask & susceptible_mask, "dynamic.vaccinatedDays"] = 0
    return result