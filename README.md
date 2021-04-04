# Redox_stability

## Description
This photochemistry model is originally developed by Chaffin+2017.
Original code can be found from the link below:
https://github.com/planetarymike/chaffin_natgeo_mars_photochemistry

The code is used in Koyama+21(ApJ, inpress)

This code is modified by Koyama to investigate the self-regulation of H & O
escape from atmospheres of early Mars.
Calculating self-regulation (redox re-equilibration) timescale for various
atmospheric conditions of early Mars.

Parameters:
- Atmospheric surface pressure (pCO2),
- surface temperature,
- O escape rate

## Requirements
- Julia 0.6.4 (sorry for my laziness not updated)
- Python3 (only for plot)

Packages
- PyPlot
- HDF5
- JLD

## Installation

Cloning in any directory with the command below
```
$ git clone git@github.com:coffiego/redox_stability.git
```

## Usage
Please change the directory path to yours in the code first.
#### To get a converged atmosphere

a) Set parameters in photochemistry_auto.jl

Parameters:
- pCO2
- surface temperature
- oxygen escape
- deposition velocity
- whether include N2 and Ar

b) run the code

```
$ julia photochemistry_auto.jl
```

#### To calc time response
Converges state (initial condition) is needed before calc time response.

a) Set parameters in photochemistry_auto_run.jl \
you can change some parameters from the initial state.


b) run the code
```
$ julia run_auto.jl
```
Output data:
Density profile, reaction rate, flux, escape, depositions, the ratio of H/O loss over time. \
Summary of time responses is saved in summary_reg_Tmod.csv


#### plot results
You just need to run python codes of plot_stat.py and make figures with appropriate functions
It is easier to use interactive mode in some IDEs like Spyder, Jupyter,,,
