Next Steps for SIRV Inspired Model Analysis
----
The next steps are to create an 'analyis.py' file within the Analysis folder that can flexibly run experiments to examine both intrapopulation and interpopulation differences in disease dynamics. The idea is to run various simulation runs and plot the SIRV values, and calculate the effective reproduction rate and reproduction rate during the free growth period. 

High level Idea: 
Generate an initialized population, simulate their progression, and return and plot their results. The key here is that we would like to run experiments over different possibly initialization than the ones already set in the config. Do not change the config but try to find a way to run different populations in a programmatic way.

Attributes we would like to vary:
Distribution of Age Categories in population
Distribution of Vaccine Acceptance in population
Distributions of Comorbidities in Population
Variance of an individual's movement (currently sig =4 in generator)
The density of the population (rho)

For each of these we would like to be able to run simulations with different values and see how the free growth r0 and effective re changes. The `analysis.py` file should fit in seamlessly with the existing codebase structure. 