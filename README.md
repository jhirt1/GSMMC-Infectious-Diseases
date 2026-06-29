[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20659487.svg)](https://doi.org/10.5281/zenodo.20659487)

Objected Oriented Software for kinetic Monte Carlo (kMC) Simulation for Modeling Coronavirus Spread
-----
This file documents the planned high-level structure of the Objected Oriented Programming implementation of a kinetic Monte Carlo model for modeling the spread of coronavirus. This portion of the project builds on existing SIR, SIRV, and other kinetic Monte Carlo Models, with the ultimate goal of exploring how a flexible formulation of state allows us to explore tiered levels of complexity. 

# Code Structure

## 'main.py'
The `main.py` file will be the main analysis file, which consist of three main modules which form the backbone of the simulation and analysis. The code within `main.py` should not be edited for the purpose of new simulation runs or new parameters, and should instead take any simulation parameters as inputs to a standardized loop.

* Generation
* Simulation
* Analysis



### Generation
The generation step is the outer loop which generates a population of N individuals from a certain distribution of relevant parameters. These parameters are split between both static attributes (does not change over time) and dynamic (can change over the course of a simulation). The generation step will then call an inner simulation loop, based on a given configuration defined by a JSON schema. 


#### Functions Called
**Generator**
This function takes in a configuration schema, and outputs an initial population for the simulation. 
1. Inputs - random seed (for reproducibility), config and schema for for population parameters.
2. Outputs - Initial population for which will be simulated for $N_{sim}$ epochs. The produced dataframes are the 'Person' (static attribtues) dataframe, and the 'State' dataframe (dynamic attributes), bijectively joined by a unique id. 

#### Person - Simulation Schema - created by generator

Table 1:
| Field Name | Definition |
|---|---|
| outputLocation | Directory where the generated simulation data will be stored. |
| directory | Directory where the generator will look for files, if needed. |
| nSim | Number of iterations for the simulation to run. |
| tSpan | Time span for each epoch. |
| populationSize | Total number of individuals in the simulation. |
| initialInfected | Initial number of infected individuals at the start of the simulation. |
| seed | Random seed used for reproducibility. |
| rho | Population density. |
| rcd | Number of days required to recover from infection. |
| ved | Vaccine efficacy delay in days. |
| vid | Vaccine efficacy with infection immunity delay in days. |
| rimmd | Number of days of immunity after recovery. |
| generatorDistribution | Distribution strategy used for generating synthetic individual attributes. |
| syntheticPopulation | Explicit sampling assumptions used to generate static, schema-valid synthetic individuals. Dynamic disease, vaccination, recovery, and location states are computed by the simulation model rather than sampled here. |
| syntheticPopulation.static | Sampling assumptions for static individual attributes. |
| syntheticPopulation.static.age | Sampling weights for individual age groups. |
| syntheticPopulation.static.comorbidity | Probability that an individual has one or more comorbidities. |
| syntheticPopulation.static.socialActivity | Sampling weights for level of social activity. |
| syntheticPopulation.static.geography | Sampling weights for geographic setting. |
| syntheticPopulation.static.mobility | Sampling weights for mobility status. |
| syntheticPopulation.static.vaccineAcceptance | Probability that an individual is willing to accept vaccination. |
| infectionRiskMultipliers | Relative risk multipliers used to modify the baseline probability of infection according to static individual attributes. A value of 1.0 represents the reference condition, 0.75 represents reduced relative risk, and 1.25 represents increased relative risk. |
| infectionRiskMultipliers.age | Infection risk multipliers by age group. |
| infectionRiskMultipliers.comorbidity | Infection risk multipliers based on whether the individual has one or more comorbidities. |
| infectionRiskMultipliers.socialActivity | Infection risk multipliers by level of social activity. |
| infectionRiskMultipliers.geography | Infection risk multipliers by geographic setting. |
| infectionRiskMultipliers.mobility | Infection risk multipliers by mobility status. |
| infectionRiskMultipliers.vaccineAcceptance | Infection risk multipliers based on whether the individual is willing to accept vaccination. |

### Simulation
The simulation phase consists of a nested loop, corresponding to $N_{sim}$ iterations of epochs of length $T_{span}$.

Each time step of an epoch consists of four actions: 

For $T < T_{span}$

1. Generate update/space: update any dynamic stored days, and move each person according to their mobility distribution.
2. Compute Probabilities: For each possible transition, compute the probabilities associated with each transition reaction.
3. Initiate Transition Reactions: Calling the reactions based on associated probabilities with reactions (if multiple options, take higher)
4. $T \mathrel{{+}{=}} 1$ and repeat from step 1.

Each simulation from $T=0 \to T= T_{span}$ corresponds to one epoch, and the results from from all $N_{sim}$ are then averaged and used for resulting analysis. 

#### Functions Called
**Reactions**
This function updates the state of the population. It does not compute the probability of each reaction, but simply updates it. 
1. Inputs - Dataframe of population 'State' at time $t$
2. Outputs - Dataframe of population 'State' at time t+1

$\textit{Possible Reactions}$
- $S \to I$: Susceptible to Infected.
- $I \to R$: Infected to Recovered.
- $R \to S$: Recovered to Susceptible.
- $S \to V$: Susceptible to Vaccinated.
- $V \to S$: Vaccinated to Susceptible. 

**Probabilities/Transition Conditions**
The transition from one state to another through a reaction channel is governed by either deterministic or probabilistic propensities, which are computed based on the current state as well as probability distributions parametrized by both dynamic and static attributes.

In the current SIVR configuration:

$$S \to I$$

Transition from susceptible to infected is governed by the total infection probability, or formally the probability of infection $P_{i}$ in a given time unit for susceptible indiviual i:

$$ P_{i} = int\{\sum_{j \in \mathcal{I}} f_{ij}P_{ij}(r_{ij}) + \xi\},$$
$$ P_{ij}(r_{ij}) := \frac{1}{\sqrt{2\pi\sigma_{r}^{2}}}\exp\{-\frac{r_{ij}^{2}}{2\sigma_{r}^{2}}\}$$

where 

- $r_{ij}:=$ the euclidean distance between individual $i$ and individual $j$
- $\sigma_{r}^{2} := $ The standard deviation for the infection probability distribution, often standardized to 2.4m
- $f_{i}:=$ A risk multiplier term taking into account the static attributes such as an individual's social activity, age, geography, mobility, comorbidities, etc. 
- $\xi :=$ A uniform random variable on $[0,1]$

The risk multipliers are given by:

infectionRiskMultipliers:
  age:
    "0-17": 0.75
    "18-49": 1.0
    "50-64": 1.0
    "65+": 1.25

  comorbidity:
    "true": 1.25
    "false": 1.0

  socialActivity:
    high: 1.25
    medium: 1.0
    low: 0.75

  geography:
    urban: 1.25
    rural: 1.0

  mobility:
    independent: 1.0
    assisted: 1.0
    immobile: 0.75

  vaccineAcceptance:
    "true": 0.75
    "false": 1.0
$$ I \to R $$
The transition from infection to recovery is in the first iteration governed by the length of time after which an infected individual recovers. This time length is governed by the variable `rcd`.

$$ R \to S $$
The transition from recovered to susceptible is governed by the length of time for which they have been recovered, where individuals become susceptible again after `sd` days.

$$ S \to V $$
Susceptible individuals are transition to vaccinated based on their demographic and behavioral multipliers, specifically:
* `age`
* `comorbidity`
* `socialActivity`
* `geography`
* `mobility`
* `vaccineAcceptance`

The multipliers for each factor is taken as an input dictionary `multipliers`.

$$ V \to S $$
The transition from vaccinated to susceptible is governed by the length of time for which they have been vaccinated, where individuals become susceptible again after `vd` days.




#### Simulation Schema
Table 2.

| Field Name | Definition |
|---|---|
| static | Static demographic, geographic, and behavioral characteristics. |
| static.guid | Unique identifier for the individual. |
| static.age | Age range of the individual. |
| static.comorbidity | Indicates whether the individual has one or more comorbidities. |
| static.socialActivity | General level of social activity. |
| static.geography | Geographic setting of the individual. |
| static.mobility | Mobility status of the individual. |
| static.vaccineAcceptance | Indicates whether the individual is willing to accept vaccination. |
| dynamic | Dynamic health, vaccination, infection, recovery, and location variables that may change over time. |
| dynamic.vaccineStatus | Indicates whether the individual is currently vaccinated. |
| dynamic.proactiveVaccine | Indicates whether the individual proactively received a vaccination. |
| dynamic.numberOfInfections | Total number of infections experienced by the individual. |
| dynamic.sirvStatus | Current SIRV status of the individual: susceptible, infected, recovered, or vaccinated. |
| dynamic.infectedDays | Number of days the individual has been infected. |
| dynamic.vaccinatedDays | Number of days since the individual was vaccinated. |
| dynamic.recoveredDays | Number of days since the individual recovered from infection. |
| dynamic.currentLocation | Current x and y coordinate location of the individual. |
| dynamic.currentLocation.xcor | Current x-coordinate location of the individual. |
| dynamic.currentLocation.ycor | Current y-coordinate location of the individual. |

#### Simulation Output
!ToDo

---

## Citation

If you use this software in academic work, please cite:

Babou, T., Capece, J., Gideon, U., Hirt, J., Le, J., Tamanna, T., Yang, H., Youn, S., & Yu, X. (2026). Software for kinetic Monte Carlo (kMC) Simulation for Modeling Coronavirus Spread (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.20659487

### BibTeX

```bibtex
@misc{babou_2026_20659487,
  author       = {Babou, Talla and
                  Capece, Julianna and
                  Gideon, Uzochi and
                  Hirt, Juliana and
                  Le, Justin and
                  Tamanna, Tamanna and
                  Yang, Haoru and
                  Youn, SangEun and
                  Yu, Xinyue},
  title        = {Software for kinetic Monte Carlo (kMC) Simulation
                   for Modeling Coronavirus Spread
                  },
  month        = jun,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.2},
  doi          = {10.5281/zenodo.20659487},
  url          = {https://doi.org/10.5281/zenodo.20659487},
}
```

---

## License

MIT License

See the `LICENSE` file for full license text.
