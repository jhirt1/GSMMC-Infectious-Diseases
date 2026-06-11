import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

def calculate_re(df: pd.DataFrame,
                 time_col: str = 't',
                 status_cols: list = ['S', 'I', 'R', 'V'],
                 id_col: str = 'static.guid',
                 rcd: int = 14,
                 ved: int = 60) -> pd.DataFrame:
    '''
    Time-varying effective reproduction number R_e(t) via the renewal (Cori)
    equation. R_e(t) = incidence(t) / sum_tau incidence(t-tau) * w(tau).

    Assumes I is active prevalence; deterministic rcd-day infectious period,
    contagious throughout (flat generation interval); vaccination drawn only
    from S. Waning (V->S after ved days) needs no correction -- it is netted
    out in the observed V series, so ved is accepted only for signature parity.
    Note: R_e is undefined (NaN) for the first rcd days (estimator burn-in).
    '''
    tot = df.groupby(time_col)[status_cols].sum().sort_index()
    S, V = tot['S'], tot['V']

    # new infections/day = S outflow not explained by the net change in V
    inc = (-S.diff() - V.diff()).clip(lower=0).fillna(0)

    w = np.full(rcd, 1.0 / rcd)            # flat GI over the rcd-day window
    iv = inc.to_numpy()
    re = np.full(len(iv), np.nan)
    for k in range(rcd, len(iv)):
        foi = iv[k - rcd:k][::-1] @ w      # with flat w, this is mean of last rcd days
        if foi > 0:
            re[k] = iv[k] / foi

    return pd.DataFrame({time_col: tot.index, 'incidence': iv, 'Re': re})


def calculate_r0(df: pd.DataFrame,
                 time_col: str = 't',
                 status_cols: list = ['S', 'I', 'R', 'V'],
                 id_col: str = 'static.guid',
                 rcd: int = 14,
                 fit_window: tuple = (1, 14)) -> float:
    '''
    R0 from the free-growth (fully susceptible) phase. Fit the exponential
    growth rate r of active cases over the unrestricted window, then invert
    Euler-Lotka for a flat rcd-day infectious period:
        R0 = rcd * r / (1 - exp(-rcd * r)).
    Set fit_window to your heating phase (before restrictions / S depletion).
    '''
    tot = df.groupby(time_col)[status_cols].sum().sort_index()
    I = tot['I']
    lo, hi = fit_window
    seg = I[(I.index >= lo) & (I.index <= hi) & (I > 0)]
    r = np.polyfit(seg.index.to_numpy(), np.log(seg.to_numpy()), 1)[0]
    return rcd * r / (1.0 - np.exp(-rcd * r))