import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
    # Ensure the DataFrame is sorted by time
    df_sorted = df.copy().sort_values(by=time_col)

    # Create a pivot table to count the number of individuals in each SIRV status at each time step and facet
    pivot_df = df_sorted.pivot_table(index=[time_col, facet_col], columns=status_col, aggfunc='size', fill_value=0).reset_index()
    facets, n_facets = df_sorted[facet_col].unique(), df_sorted[facet_col].nunique()
    # Plotting with seaborn
    for i, facet in enumerate(facets):
        facet_df = df_sorted[df_sorted[facet_col] == facet]
        plot_sirv(facet_df, title=f'{title} - {facet}', id_col=id_col, status_col=status_col, time_col=time_col)