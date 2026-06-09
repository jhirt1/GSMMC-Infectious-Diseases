# GSMMC Day 1

## Code Review
* KMC_SIR_Driver
    * Each team member will need a unique seed to ensure local repeatability
    * Instead of using a random generator for the distributions, we will need to adopt a standard distribution.
    * ```tspan``` is an abstract variable, which will need to be interpreted by the team.

* SusceptibleToInfected
    * Compute the distances, compute the probability, profit.

* Questions
    * Are the recovery days checked at the beginning or end of an epoch?
        * Technically it can be done at the beginning or end. As of right now it's at the beginning.
    * Does this model follow the basic assumptions of probability?
        * Not right now!!! 

