import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_sirv_onehot(df: pd.DataFrame,
            title: str = 'SIRV Status Over Time',
            id_col: str = 'static.guid',
            run_col: str = 'run',
            time_col: str = 't') -> None:
    '''
    Alternate plotting code for one-hot encoded SIRV status.
    '''
    df_sorted = df.copy().sort_values(by=time_col)

    
    status_cols = ['S', 'I', 'R', 'V']
    # Use the agg function to count mean number of each status per time point
    agg_df = df_sorted.groupby(time_col)[status_cols].sum().reset_index()
    plt.figure(figsize=(12, 6))
    for status in status_cols:
        plt.plot(agg_df[time_col], agg_df[status], label=status)
    plt.xlabel('Time')
    plt.ylabel('Count')
    plt.title(title)
    plt.legend()
    plt.show()

def plot_sirv_onehot_by_facet(df: pd.DataFrame,
                            facet_col: str,
                            id_col: str = 'static.guid',
                            time_col: str = 't') -> None:
    '''
    Plots the one-hot encoded SIRV status of individuals over time, faceted by a specified column.
    '''
    df_sorted = df.copy().sort_values(by=time_col)
    facets, n_facets = df_sorted[facet_col].unique(), df_sorted[facet_col].nunique()
    for i, facet in enumerate(facets):
        facet_df = df_sorted[df_sorted[facet_col] == facet]
        plot_sirv_onehot(facet_df, title= f'SIRV by Facet - {facet}', id_col=id_col, time_col=time_col)

def plot_sirv(df: pd.DataFrame,
            title: str = 'SIRV Status Over Time',
            id_col: str = 'static.guid',
            status_col: str = 'dynamic.sirvStatus',
            time_col: str = 't') -> None:
    '''
    Plots the SIRV status of individuals over time.
    '''
    # Ensure the DataFrame is sorted by time
    df_sorted = df.copy().sort_values(by=time_col)

    # Create a pivot table to count the number of individuals in each SIRV status at each time step
    pivot_df = df_sorted.pivot_table(index=time_col, columns=status_col, aggfunc='size', fill_value=0)

    # Plotting
    plt.figure(figsize=(12, 6))
    for status in df[status_col].unique():
        if status in pivot_df.columns:
            plt.plot(pivot_df.index, pivot_df[status], label=status)
    plt.title(title)
    plt.xlabel('Time')
    plt.ylabel('Number of Individuals')
    plt.legend()
    plt.show()

def plot_sirv_by_facet(df: pd.DataFrame,
            facet_col: str,
            title: str = 'SIRV Status Over Time by Facet',
            id_col: str = 'static.guid',
            status_col: str = 'dynamic.sirvStatus',
            time_col: str = 't') -> None:
    '''
    Plots the SIRV status of individuals over time, faceted by a specified column.
    '''
    df_copy = df.copy()
    facets, n_facets = df_copy[facet_col].unique(), df_copy[facet_col].nunique()
    # Plotting
    for i, facet in enumerate(facets):
        facet_df = df_copy[df_copy[facet_col] == facet]
        plot_sirv(facet_df, title=f'{title} - {facet}', id_col=id_col, status_col=status_col, time_col=time_col)