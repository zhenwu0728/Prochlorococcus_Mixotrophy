# README

[![DOI](https://zenodo.org/badge/438721739.svg)](https://zenodo.org/badge/latestdoi/438721739)


This repository contains the code and scripts for paper "Significant organic carbon acquisition by Prochlorococcus"

## Growth rate estimation
This folder contains the code to evaluate the potential growth rates at the Equatorial Pacific (EqPac, 0oN, 140oW) and North Pacific Subtropical Gyre (HOT, 22o45’N, 158oW; Station ALOHA). The Jupyter Notebook `figure2.ipynb` runs the scripts and generate Figure 2 in the paper. `Julia v1.6.5` is required to run this notebook. The `Julia` packages used in this notebook include `CSV.jl`, `DataFrames.jl`, `Polynomials.jl`, `Printf.jl`, `Parameters.jl`, `Statistics.jl`. We use `CairoMakie.jl`, and `LaTeXStrings.jl` for plotting.

## Individual-based simulation
This folder contains the script and parameters to run the individual-based model in the paper. `Julia v1.6.5` is required to run the model. We use `PlanktonIndividuals.jl v0.1.9` to simulate phytoplankton cells. We use `Oceananigans.jl v0.55.0` to generate the flow fields. Extra `Julia` packages used in this script include `Random.jl`, `Statistics.jl`, `Printf.jl`, `YAML.jl`, `Serialization.jl`, `SeawaterPolynomials.jl`.

The simulation is run for 360 simulated days with a time step of 1 minute. The whole simulation takes about 1 day to finish. The model outputs are binary files, we suggest to use `Serialization.jl` to read these output files.
