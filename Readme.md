# README

[![DOI](https://zenodo.org/badge/438721739.svg)](https://zenodo.org/badge/latestdoi/438721739)


This repository contains the code and scripts for paper "Significant organic carbon acquisition by *Prochlorococcus*"

## Growth rate estimation
This folder contains the code to evaluate the potential growth rates at the Equatorial Pacific (EqPac, 0N, 140W) and North Pacific Subtropical Gyre (HOT, 22.75N, 158W; Station ALOHA). The input files for this evaluation are `hourly_PAR.csv` (hourly light intensity) and `pro_mixo_input.csv` (nutrient concentrations, cell densities and growth rates of *Prochlorococcus*). The growth rates of *Prochlorococcus* were from Vaulot et al. and Liu et al. The concurrent photon fluxes and nutrient concentrations were available from an extensive biogeochemical survey (JGOFS EqPac) and time-series station (HOT, cruise #55).

The Jupyter Notebook `figure2.ipynb` runs the scripts and generate Figure 2 in the paper. `Julia v1.6.5` is required to run this notebook. The `Julia` packages used in this notebook include `CSV.jl`, `DataFrames.jl`, `Polynomials.jl`, `Printf.jl`, `Parameters.jl`, `Statistics.jl`. We use `CairoMakie.jl`, and `LaTeXStrings.jl` for plotting.

### References

Vaulot, D., Marie, D., Olson, R. J. & Chisholm, S. W. Growth of Prochlorococcus, a photosynthetic prokaryote, in the equatorial Pacific Ocean. Science (1979) 268, 1480–1482 (1995).

Liu, H., Nolla, H. & Campbell, L. Prochlorococcus growth rate and contribution to primary production in the equatorial and subtropical North Pacific Ocean. Aquatic Microbial Ecology 12, 39–47 (1997).

Murray, J., Leinen, M., Feely, R., Toggweiler, R. & Wanninkhof, R. EqPac: A Process Study in the Central Equatorial Pacific. Oceanography 5, 134–142 (1992).

Karl, D. M. & Church, M. J. Microbial oceanography and the Hawaii Ocean Time-series programme. Nature Reviews Microbiology 12, 699–713 (2014)


## Individual-based simulation
This folder contains the script and parameters to run the individual-based model in the paper. `Julia v1.6.5` is required to run the model. We use [`PlanktonIndividuals.jl v0.1.9`](https://github.com/JuliaOcean/PlanktonIndividuals.jl) to simulate phytoplankton cells. We use `Oceananigans.jl v0.55.0` to generate the flow fields. Extra `Julia` packages used in this script include `Random.jl`, `Statistics.jl`, `Printf.jl`, `YAML.jl`, `Serialization.jl`, `SeawaterPolynomials.jl`.

The simulation is run for 360 simulated days with a time step of 1 minute. The whole simulation takes about 1 day to finish. The initial conditions of the simulation are binary files located in `Input` folder. The model outputs are also binary files, we suggest to use `Serialization.jl` to read these output files.
