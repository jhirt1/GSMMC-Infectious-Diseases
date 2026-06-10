import pandas as pd
import numpy as np

def infected_to_recovered(df: pd.DataFrame, rcd: int) -> pd.DataFrame:
    '''
    Converts an infected individual to recovered after rcd days of infection. 
    The individual does not join the susceptible pool
    '''
    result = df.copy()
    
    # Find infected individuals who have reached recovery day
    infected_mask = result['sirvStatus'] == 'I'
    recovered_mask = result['infectedDays'] == rcd
    
    # Update status to 'recovered' for those who've completed rcd days
    result.loc[infected_mask & recovered_mask, 'sirvStatus'] = 'R'
    
    return result

def infection_probability(df: pd.DataFrame, sig2: float) -> np.ndarray:
    '''
    Compute infection probability for each susceptible individual based on proximity and infectivity of infected individuals.
    Returns an array of probabilities aligned to the full DataFrame index.
    '''
    tsig2 = 2 * sig2

    susceptible_mask = df['sirvStatus'] == 'S'
    infected_mask = df['sirvStatus'] == 'I'

    susceptible_idx = df.index[susceptible_mask]
    infected_idx = df.index[infected_mask]

    pi_full = np.zeros(len(df))

    if len(susceptible_idx) == 0 or len(infected_idx) == 0:
        return pi_full

    xi = df.loc[susceptible_idx, 'x'].values
    yi = df.loc[susceptible_idx, 'y'].values

    xj = df.loc[infected_idx, 'x'].values
    yj = df.loc[infected_idx, 'y'].values
    aRM = df.loc[infected_idx, 'ageRiskMultiplier'].values
    cRM = df.loc[infected_idx, 'comorbidityRiskMultiplier'].values
    saRM = df.loc[infected_idx, 'socialActivityRiskMultiplier'].values
    gRM = df.loc[infected_idx, 'geographyRiskMultiplier'].values
    mRM = df.loc[infected_idx, 'mobilityRiskMultiplier'].values
    vaRM = df.loc[infected_idx, 'vaccineAcceptanceRiskMultiplier'].values

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

    susceptible_mask = result['sirvStatus'] == 'S'
    susceptible_positions = np.where(susceptible_mask.values)[0]

    if len(susceptible_positions) == 0:
        return result

    # Stochastic threshold: infect if floor(pi + U[0,1)) >= 1
    pi_susceptible = pi[susceptible_positions]
    newly_infected = np.floor(pi_susceptible + np.random.rand(len(susceptible_positions))).astype(bool)

    newly_infected_idx = result.index[susceptible_positions[newly_infected]]
    result.loc[newly_infected_idx, 'sirvStatus'] = 'I'

    return result

def recovered_to_susceptible(df: pd.DataFrame, sd: int) -> pd.DataFrame:
    '''
    Converts a recovered individual to susceptible after sd days. 
    '''
    result = df.copy()

    # Find recovered individuals who have reached susceptibility day
    recovered_mask = result['sirvStatus'] == 'R'
    susceptible_mask = result['recoveredDays'] == sd

    # Update status to 'susceptible' for those who've completed sd days
    result.loc[recovered_mask & susceptible_mask, 'sirvStatus'] = 'S'
    return result

def susceptible_to_vaccinated(df: pd.DataFrame, multipliers: dict, base_prob: float = 0.005) -> pd.DataFrame:
    '''
    Converts susceptible individuals to vaccinated based on demographic and behavioural multipliers.
    '''
    result = df.copy()

    susceptible_mask = result['sirvStatus'] == 'S'
    susceptible_positions = np.where(susceptible_mask.values)[0]

    if len(susceptible_positions) == 0:
        return result

    sus = result.iloc[susceptible_positions]

    prob = (
        base_prob
        * sus['age'].map(multipliers['age'])
        * sus['comorbidity'].map(multipliers['comorbidity'])
        * sus['socialActivity'].map(multipliers['socialActivity'])
        * sus['geography'].map(multipliers['geography'])
        * sus['mobility'].map(multipliers['mobility'])
        * sus['vaccineAcceptance'].map(multipliers['vaccineAcceptance'])
    ).clip(0, 1).values

    newly_vaccinated = np.random.rand(len(susceptible_positions)) < prob

    newly_vaccinated_idx = result.index[susceptible_positions[newly_vaccinated]]
    result.loc[newly_vaccinated_idx, 'sirvStatus'] = 'V'

    return result

def vaccinated_to_susceptible(df: pd.DataFrame, vd: int) -> pd.DataFrame:
    '''
    Converts a vaccinated individual to susceptible after vd days of vaccination. 
    '''
    result = df.copy()

    # Find vaccinated individuals who have reached susceptibility day
    vaccinated_mask = result['sirvStatus'] == 'V'
    susceptible_mask = result['vaccinatedDays'] == vd

    # Update status to 'susceptible' for those who've completed sd days
    result.loc[vaccinated_mask & susceptible_mask, 'sirvStatus'] = 'S'
    return result
